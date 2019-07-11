/* EnvironentSensing

   For KPIC FIU and HCST-R.

   Feather M0 Adalogger w/ Ethernet FeatherWing
   connected via I2C to TCA9548A I2C multiplexer
   Each I2C channel on the mux has one BME680 and one LSM9DS1.
   
   Modified July 2, 2019
   by Grady Morrissey
   
   Modified August 23, 2018
   by Milan Roberson
*/

#include "Wire.h"
#include "Adafruit_BME680.h"
#include "Adafruit_LSM9DS1.h"
#include "SPI.h"
#include "Ethernet2.h"
#include "EthernetUdp2.h"
#include "SdFat.h"
#include "TimeLib.h"
#include "avr/dtostrf.h"
#include "math.h"

#define _TASK_PRIORITY
#define _TASK_WDT_IDS
#include "TaskScheduler.h"

#include <Timezone.h> 

#define MUXADDR 0x74

#define NUMPORTS 8

#define CADENCE_SECONDS 1

//#define DEBUG

#ifdef DEBUG
#define Sprint(a) (Serial.print(a))
#define Sprintln(a) (Serial.println(a))
#define Sbegin(a) {Serial.begin(115200); while (!Serial) delay(1000);}
#else
#define Sprint(a)
#define Sprintln(a)
#define Sbegin(a)
#endif

Adafruit_BME680 bme[NUMPORTS];
Adafruit_LSM9DS1 lsm[NUMPORTS];

// Web server globals.
// mac addresss is unique to every Ethernet FeatherWing.
//GRADY NOTE
uint8_t mac[] = {0x98, 0x76, 0xB6, 0x10, 0xB4, 0x1A};
EthernetServer server(80);


// NTP globals
// NTP Servers:
//GRADY NOTE
IPAddress timeServer(132, 163, 96, 1); 
//IPAddress timeServer(132, 163, 97, 1); // time-a-wwv.nist.gov
//IPAddress timeServer(132, 163, 97, 2); // time-b-www.nist.gov
// IPAddress timeServer(132, 163, 97, 3); // time-c-www.nist.gov

//GRADY NOTE
// const int timeZone = -8;  // Pacific Standard Time (USA)
// const int timeZone = -7;  // Pacific Daylight Time (USA)
// const int timeZone = -10; // Hawaii Standard Time (USA)
// const int timeZone = 0; // UTC


// US Pacific Time Zone (Las Vegas, Los Angeles)
TimeChangeRule usPDT = {"PDT", Second, Sun, Mar, 2, -420};
TimeChangeRule usPST = {"PST", First, Sun, Nov, 2, -480};
Timezone local(usPDT, usPST);





EthernetUDP Udp;
const unsigned int localPort = 8888;


// Data is held here
bool isSetUp[NUMPORTS];
float temps[NUMPORTS];
float hums[NUMPORTS];
float pressure[NUMPORTS];
float accx[NUMPORTS];
float accy[NUMPORTS];
float accz[NUMPORTS];
float gyrox[NUMPORTS];
float gyroy[NUMPORTS];
float gyroz[NUMPORTS];
float magx[NUMPORTS];
float magy[NUMPORTS];
float magz[NUMPORTS];

// SD card stuff
#define cardSelect 4

char rootFileName[] = "index.htm";

SdFat card;
SdFile file;

void handleRequest();
void takeMeasurements();

Scheduler r, hpr;
Task serve(TASK_SECOND / 2, TASK_FOREVER, &handleRequest, &r);
Task sense(TASK_SECOND * CADENCE_SECONDS, TASK_FOREVER, &takeMeasurements, &hpr);

void tcaselect(uint8_t i) {
  if (i > 7) return;

  Wire.beginTransmission(MUXADDR);
  Wire.write(1 << i);
  Wire.endTransmission();
}

/* NTP code */
const int NTP_PACKET_SIZE = 48; // NTP time is in the first 48 bytes of message
byte packetBuffer[NTP_PACKET_SIZE]; //buffer to hold incoming & outgoing packets

time_t getNtpTime() {
  while (Udp.parsePacket() > 0) ; // discard any previously received packets
  Sprintln("Transmit NTP Request");
  sendNTPpacket(timeServer);
  uint32_t beginWait = millis();
  while (millis() - beginWait < 1500) {
    int size = Udp.parsePacket();
    if (size >= NTP_PACKET_SIZE) {
      Sprintln("Receive NTP Response");
      Udp.read(packetBuffer, NTP_PACKET_SIZE);  // read packet into the buffer
      unsigned long secsSince1900;
      // convert four bytes starting at location 40 to a long integer
      secsSince1900 =  (unsigned long)packetBuffer[40] << 24;
      secsSince1900 |= (unsigned long)packetBuffer[41] << 16;
      secsSince1900 |= (unsigned long)packetBuffer[42] << 8;
      secsSince1900 |= (unsigned long)packetBuffer[43];
      return secsSince1900 - 2208988800UL;
    }
  }
  Sprintln("No NTP Response :-(");
  return 0; // return 0 if unable to get the time
}

// send an NTP request to the time server at the given address
void sendNTPpacket(IPAddress &address) {
  // set all bytes in the buffer to 0
  memset(packetBuffer, 0, NTP_PACKET_SIZE);
  // Initialize values needed to form NTP request
  // (see URL above for details on the packets)
  packetBuffer[0] = 0b11100011;   // LI, Version, Mode
  packetBuffer[1] = 0;     // Stratum, or type of clock
  packetBuffer[2] = 6;     // Polling Interval
  packetBuffer[3] = 0xEC;  // Peer Clock Precision
  // 8 bytes of zero for Root Delay & Root Dispersion
  packetBuffer[12]  = 49;
  packetBuffer[13]  = 0x4E;
  packetBuffer[14]  = 49;
  packetBuffer[15]  = 52;
  // all NTP fields have been given values, now
  // you can send a packet requesting a timestamp:
  Udp.beginPacket(address, 123); //NTP requests are to port 123
  Udp.write(packetBuffer, NTP_PACKET_SIZE);
  if (!Udp.endPacket()) {
    Sprintln("Error sending UDP packet");
  }
}

void printHeader(Print &f) {
  char text[141];
  digitalWrite(8, HIGH); // turns LED on
  f.write("# Time");
  for (uint8_t i = 0; i < NUMPORTS; i++) {
    snprintf(text, 140, ", Temp %d, Humidity %d, Pressure %d, Accel X %d, Accel Y %d, Accel Z %d, Gyro X %d, Gyro Y %d, Gyro Z %d, Mag X %d, Mag Y %d, Mag Z %d", i, i, i, i, i, i, i, i, i, i, i, i);
    f.write(text);
  }
  f.write("\n# Time is UTC, Temp is *C, Humidity is %, Pressure is hPa or millibar, Accel is m/s^2, Gyro is degrees/s, Mag is gauss\n");
  digitalWrite(8, LOW);
}

void getFileName(char *filename, time_t t) {
  snprintf(filename, 13, "%04d%02d%02d.csv", year(t), month(t), day(t));
}

void writeData(time_t t, Print &f) {
  char data[11];

  digitalWrite(8, HIGH); // turns LED on

  snprintf(data, 10, "%02d:%02d:%02d", hour(t), minute(t), second(t));
  f.write(data);
  for (uint8_t i = 0; i < NUMPORTS; i++) {
    f.write(",");
    dtostrf(temps[i], 6, 2, data);
    f.write(data);
    f.write(",");
    dtostrf(hums[i], 5, 2, data);
    f.write(data);
    f.write(",");
    dtostrf(pressure[i], 6, 2, data);
    f.write(data);
    f.write(",");
    dtostrf(accx[i], 6, 3, data);
    f.write(data);
    f.write(",");
    dtostrf(accy[i], 6, 3, data);
    f.write(data);
    f.write(",");
    dtostrf(accz[i], 6, 3, data);
    f.write(data);
    f.write(",");
    dtostrf(gyrox[i], 6, 2, data);
    f.write(data);
    f.write(",");
    dtostrf(gyroy[i], 6, 2, data);
    f.write(data);
    f.write(",");
    dtostrf(gyroz[i], 6, 2, data);
    f.write(data);
    f.write(",");
    dtostrf(magx[i], 6, 3, data);
    f.write(data);
    f.write(",");
    dtostrf(magy[i], 6, 3, data);
    f.write(data);
    f.write(",");
    dtostrf(magz[i], 6, 3, data);
    f.write(data);
  }
  f.write("\n");
  digitalWrite(8, LOW);
}

void setup() {

  // for debugging
  Sbegin(115200);

  // Set up communication with Mux
  Wire.begin();
  Sprintln("I2C mux scanner ready");

  // Check that sensors are connected
  for (uint8_t t = 0; t < NUMPORTS; t++) {
    tcaselect(t);
    Sprint("mux port #"); Sprintln(t);
    if (bme[t].begin()) {
      bme[t].setGasHeater(0, 0); // disable gas reading

      lsm[t].begin();
      isSetUp[t] = true;
      Sprint("Setting up sensors "); Sprintln(t);
    } else {
      isSetUp[t] = false;
    }
  }
  
  // Needed for Ethernet
  pinMode(13, OUTPUT);
  pinMode(8, OUTPUT);

  // Set up SD card
  pinMode(10, OUTPUT);
  digitalWrite(10, HIGH);
  if (!card.begin(cardSelect, SPI_FULL_SPEED)) {
    Sprintln("card.init failed");
    setup();
  } else {
    card.vwd()->rewind();
    card.chdir();
    Sprintln("Files found in root:");
#ifdef DEBUG
    card.ls(LS_DATE | LS_SIZE);
    Sprintln();
    Sprintln("Files found in all dirs:");
    card.ls(LS_R);
    Sprintln();
    Sprintln("Done");
#endif
    FatFile::dateTimeCallback(sd_dateTimeCallback);

    file.open("config.txt", O_READ);
    char ip_s[25];

    file.read(ip_s, 12);
    if (ip_s[11] == '1') {
      // read in string
      file.read(ip_s, 24);
      IPAddress gateway;
      gateway.fromString((ip_s + 9));

      file.read(ip_s, 24);
      IPAddress subnet;
      subnet.fromString((ip_s + 9));

      file.read(ip_s, 20);
      IPAddress ip;
      ip_s[20] = '\0';
      ip.fromString((ip_s + 5));
      Ethernet.begin(mac, ip, gateway, subnet);
    } else {
      Ethernet.begin(mac);
    }
    
    file.close();
  }

  server.begin();
  Sprint("server is at "); Sprintln(Ethernet.localIP());

  // Set up NTP
  Udp.begin(localPort);
  setSyncProvider(getNtpTime);

  r.setHighPriorityScheduler(&hpr);
  r.enableAll(true);
}

void loop() {
  r.execute();
}

time_t global_t;
void takeMeasurements() {
  global_t = now();
  global_t = local.toLocal(global_t);
  for (uint8_t i = 0; i < NUMPORTS; i++) {
    tcaselect(i);
    if (!isSetUp[i] && bme[i].begin()) {
      bme[i].setGasHeater(0, 0); // disable gas reading
      lsm[i].begin();
    }
    isSetUp[i] = bme[i].performReading();
    if (isSetUp[i]) {
      temps[i] = bme[i].temperature;
      hums[i] = bme[i].humidity;
      pressure[i] = bme[i].pressure / 100; // convert Pa to hPa

      lsm[i].read();
      sensors_event_t a, m, g, temp;
      lsm[i].getEvent(&a, &m, &g, &temp);
      accx[i] = a.acceleration.x;
      accy[i] = a.acceleration.y;
      accz[i] = a.acceleration.z;

      magx[i] = m.magnetic.x;
      magy[i] = m.magnetic.y;
      magz[i] = m.magnetic.z;

      gyrox[i] = g.gyro.x;
      gyroy[i] = g.gyro.y;
      gyroz[i] = g.gyro.z;
    } else {
      temps[i] = NAN;
      hums[i] = NAN;
      pressure[i] = NAN;

      accx[i] = NAN;
      accy[i] = NAN;
      accz[i] = NAN;

      magx[i] = NAN;
      magy[i] = NAN;
      magz[i] = NAN;

      gyrox[i] = NAN;
      gyroy[i] = NAN;
      gyroz[i] = NAN;
    }
  }

  // append to today's file
  char filename[13];
  getFileName(filename, global_t);
  file.close();
  if (!file.open(filename, O_WRITE | O_APPEND)) {
    file.open(rootFileName, O_WRITE);
    file.seekEnd(-25);
    file.write("<tr>\n<td><a href='");
    file.write(filename);
    file.write("'>");
    file.write(filename);
    file.write("</a></td>\n<td>");
    snprintf(filename, 11, "%04d-%02d-%02d", year(global_t), month(global_t), day(global_t));
    file.write(filename);
    file.write("</td>\n</tr>\n</table>\n</body>\n</html>");
    file.close();
    getFileName(filename, global_t);

    file.open(filename, O_CREAT | O_WRITE | O_APPEND);
    printHeader(file);
  }
  writeData(global_t, file);
  file.sync();
  file.close();
}

void handleRequest() {
  Ethernet.maintain();
  EthernetClient client = server.available();
  if (client) {
    uint8_t bufindex = 0;
    const uint8_t maxbyte = 255;
    uint8_t buf[maxbyte];

    char *filename;
    char c;

    Sprintln("new client");
    while (client.connected()) {
      if (client.available()) {
        c = client.read();
        Sprint(c);
        if (c != '\n' && c != '\r') {
          buf[bufindex++] = c;
          if (bufindex >= maxbyte) {
            bufindex--;
          }
          continue;
        }
        buf[bufindex] = 0;
        filename = 0;
        Sprintln();
        Sprintln((char *)buf);
        if (strstr((char *)buf, "GET / ") != 0 || strstr((char *)buf, "HEAD / ") != 0) {
          filename = rootFileName;
        }
        if (strstr((char *)buf, "GET /") != 0 || strstr((char *)buf, "HEAD /") != 0) {
          if (strstr((char *)buf, "GET /now") != 0) {
            client.println("HTTP/1.1 200 OK");
            client.println("Content-type: text/csv");
            // client.println("Content-type: text/plain");
            client.println("Connection: close");
            client.println();
            printHeader(client);
            writeData(global_t, client);
          } else if (strstr((char *)buf, "GET /start") != 0) {
            client.println("HTTP/1.1 200 OK");
            client.println("Content-type: text/plain");
            client.println("Connection: close");
            client.println();
            sense.enable();
            client.println("Sensing started.");
          } else if (strstr((char *)buf, "GET /stop") != 0) {
            client.println("HTTP/1.1 200 OK");
            client.println("Content-type: text/plain");
            client.println("Connection: close");
            client.println();
            sense.disable();
            client.println("Sensing stopped.");
          } else if (strstr((char *)buf, "GET /delete") != 0) {
            client.println("HTTP/1.1 200 OK");
            client.println("Content-type: text/plain");
            client.println("Connection: close");
            client.println();
            deleteData();
            client.println("Data deleted.");
          } else {
            if (!filename) filename = (char *)buf + 5;
            (strstr((char *)buf, " HTTP"))[0] = 0;
            Sprintln(filename);
            // card.chdir();
            Sprintln(FreeRam());
            if (! file.open(filename, O_READ)) {
              Sprintln("404 not found");
              client.println("HTTP/1.1 404 Not Found");
              client.println("Content-type: text/html");
              client.println("Connection: close");
              client.println();
              if (strstr((char *)buf, "GET")) {
                file.open("404error.htm", O_READ);
                bufindex = 0;
                while ((bufindex = file.read(buf, maxbyte)) == maxbyte) {
                  client.write(buf, maxbyte);
                }
                file.close();
                client.write(buf, bufindex);
              }
              break;
            } // 404 not found
            client.println("HTTP/1.1 200 OK");
            Sprintln("Opened!");

            if (strstr(filename, ".htm") != 0) {
              client.println("Content-type: text/html");
            } else if (strstr(filename, ".css") != 0) {
              client.println("Content-type: text/css");
            } else if (strstr(filename, ".csv") != 0) {
              client.println("Content-type: text/csv");
            } else if (strstr(filename, ".ico") != 0) {
              client.println("Content-type: image/x-icon");
            } else {
              client.println("Content-type: text/plain");
            }
            client.println("Connection: close");
            client.println();
            if (strstr((char *)buf, "GET")) {
              bufindex = 0;
              while ((bufindex = file.read(buf, maxbyte)) == maxbyte) {
                client.write(buf, maxbyte);
              }
              file.close();
              client.write(buf, bufindex);
            }
            file.close();
          }
        } else {
          client.println("HTTP/1.1 404 Not Found");
          client.println("Content-type: text/html");
          client.println("Connection: close");
          client.println();
          file.open("404error.htm", O_READ);
          bufindex = 0;
          while ((bufindex = file.read(buf, maxbyte)) == maxbyte) {
            client.write(buf, maxbyte);
          }
          file.close();
          client.write(buf, bufindex);
        }
        break;
      }
    }
    delay(1);
    client.stop();
  }
}

void sd_dateTimeCallback(uint16_t* date, uint16_t* time) {
  time_t m_t = now();
  *date = FAT_DATE(year(m_t), month(m_t), day(m_t));
  *time = FAT_TIME(hour(m_t), minute(m_t), second(m_t));
}

void deleteData() {
  char fname[15];
  char cur_fname[15];
  uint8_t buf[255];
  int16_t index;

  SdFile file2;
  card.vwd()->rewind();
  card.chdir();

  getFileName(cur_fname, global_t);
  // while there are more files
  while (file.openNext(card.vwd(), O_READ)) {
    file.getName(fname, 15);
    if (strstr(fname, ".csv") != 0) {
      if (strstr(fname, cur_fname) == 0) {
        // if it is a .csv file and it's not the current file we're writing to, delete it.
        file.close();
        card.remove(fname);
      }
    }
    file.close();
  }
  // delete the old index.htm file
  file.open(rootFileName, O_RDWR);
  file.remove();
  file.close();
  file.open(rootFileName, O_WRITE | O_CREAT);
  file2.open("index_ar.htm", O_READ);
  // copy the archived index file to index.htm
  while ((index = file2.read(buf, 255)) == 255) {
    file.write(buf, 255);
  }
  file.write(buf, index);
  // write the current file into the index.
  file.seekEnd(-25);
  file.write("<tr>\n<td><a href='");
  file.write(cur_fname);
  file.write("'>");
  file.write(cur_fname);
  file.write("</a></td>\n<td>");
  snprintf(cur_fname, 11, "%04d-%02d-%02d", year(global_t), month(global_t), day(global_t));
  file.write(cur_fname);
  file.write("</td>\n</tr>\n</table>\n</body>\n</html>");
  file.close();
  file2.close();
}

extern "C" char *sbrk(int i);
int FreeRam () {
  char stack_dummy = 0;
  return &stack_dummy - sbrk(0);
}

// TODO: write documentation
