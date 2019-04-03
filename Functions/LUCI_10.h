//
//
//		Header for LUCI_10.dll
//		For use with LUCI-10 USB Interface.
//
//		This Header file is both valid for the 32- and 64-bit version of the DLL.
//	
//		Version		2.0		
//		Date		20.03.17
//
//
//		


#ifndef _FEMTO_LUCI10_H_
#define _FEMTO_LUCI10_H_

//	Error return codes
#define LUCI_OK				0	//	Function returns successful
#define LUCI_ERROR_INDEX	-1	//	Error, selected LUCI-10 not in list
#define LUCI_ERROR_HID		-2	//	Error, LUCI-10 doesnt respond

//	return value of GetStatusPin
#define	STATUS_PIN_LOW		0	//	probed Pin is logical 0	(TTL low  0 Volt)
#define STATUS_PIN_HIGH		1	//	probed Pin is logical 1	(TTL high 5 Volt)


//	name of the DLL to be loaded
#ifdef _WIN64
#define FEMTO_LUCI10_DLL "LUCI_10_x64"
#else
#define FEMTO_LUCI10_DLL "LUCI_10"
#endif

//	Implement the DLL export/import mechanism and allow a C-written program
//	to use this DLL
//#ifdef	FEMTO_LUCI10_EXPORTS
//#define FEMTO_LUCI10_API extern "C" __declspec(dllexport)
//#else
//#define FEMTO_LUCI10_API extern "C" __declspec(dllimport)
//#endif

//
//	list of exported functions
//
//FEMTO_LUCI10_API int EnumerateUsbDevices();
//FEMTO_LUCI10_API int LedOn(int index);
//FEMTO_LUCI10_API int LedOff(int index);
//FEMTO_LUCI10_API int ReadAdapterID(int index, int *id);
//FEMTO_LUCI10_API int WriteAdapterID(int index, int id);
//FEMTO_LUCI10_API int FirmwareUpdate(int index);
//FEMTO_LUCI10_API int WriteData(int index, int data_low, int data_high);
//FEMTO_LUCI10_API int GetStatusPin5(int index, int *status);
//FEMTO_LUCI10_API int GetStatusPin6(int index, int *status);
//FEMTO_LUCI10_API int GetStatusPin7(int index, int *status);
//FEMTO_LUCI10_API int GetProductString(int index, char *string, int size);

 #ifdef __cplusplus
 extern “C” {
 #endif

int EnumerateUsbDevices();
int LedOn(int index);
int LedOff(int index);
int ReadAdapterID(int index, int *id);
int WriteAdapterID(int index, int id);
int FirmwareUpdate(int index);
int WriteData(int index, int data_low, int data_high);
int GetStatusPin5(int index, int *status);
int GetStatusPin6(int index, int *status);
int GetStatusPin7(int index, int *status);
int GetProductString(int index, char *string, int size);


 #ifdef __cplusplus
 }
 #endif

#endif  // _FEMTO_LUCI10_H_