// Date and time functions using a DS1307 RTC connected via I2C and Wire lib
#include <Wire.h>
#include "RTClib.h"
#include <DS1307.h>
#include <Adafruit_RGBLCDShield.h>
RTC_DS1307 rtc;

Adafruit_RGBLCDShield lcd = Adafruit_RGBLCDShield();
#define GREEN 0x2
#define RED 0x1
#define YELLOW 0x3
#define GREEN 0x2
#define TEAL 0x6
#define BLUE 0x4
#define VIOLET 0x5
#define WHITE 0x7


char daysOfTheWeek[7][12] = {"Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"};
void setup () {
 
 Serial.begin(9600);
 lcd.begin(16, 2);              // start the library
  lcd.setBacklight(WHITE);
  rtc.begin();
 
 if (! rtc.begin()) {
   Serial.println("Couldn't find RTC");
   abort();
   
 }

 
 //rtc.adjust(DateTime(2021, 3, 11, 8, 15, 0));
 if (! rtc.isrunning()) {
   Serial.println("RTC is NOT running!");
   rtc.adjust(DateTime(F(__DATE__), F(__TIME__))); //Adjusts based on computer time
   // following line sets the RTC to the date & time this sketch was compiled
   // This line sets the RTC with an explicit date & time, for example to set
   // manually set yyyy, mm, dd, hr, m, s
   //rtc.adjust(DateTime(2021, 3, 11, 8, 26, 0));
   
 }else if (rtc.isrunning()) {
   Serial.println("RTC is running!");
   if(rtc.now()!=DateTime(F(__DATE__), F(__TIME__))){
    rtc.adjust(DateTime(F(__DATE__), F(__TIME__))); 
    Serial.println("Time set!");
    //rtc.adjust(DateTime(2021, 3, 11, 8, 26, 0));
   }
   
 }
}

void loop () {
   if (! rtc.isrunning()) {
   Serial.println("RTC is NOT running!");}
   if (rtc.isrunning()) {
   Serial.println("RTC is running!");}
   timeDisplay();
 DateTime now = rtc.now();
 Serial.println(now.timestamp());
 Serial.print(now.year(), DEC);
 Serial.print('/');
 Serial.print(now.month(), DEC);
 Serial.print('/');
 Serial.print(now.day(), DEC);
 Serial.print(" (");
// Serial.print(daysOfTheWeek[now.dayOfTheWeek()]);
 Serial.print(") ");
 Serial.print(now.hour(), DEC);
 Serial.print(':');
 Serial.print(now.minute(), DEC);
 Serial.print(':');
 Serial.print(now.second(), DEC);
 Serial.println();
 //Serial.print(" since midnight 1/1/1970 = ");
 //Serial.print(now.unixtime());
 //Serial.print("s = ");
 //Serial.print(now.unixtime() / 86400L);
 //Serial.println("d");
 // calculate a date which is 7 days and 30 seconds into the future
 //DateTime future (now + TimeSpan(7, 12, 30, 6));
 //Serial.print(" now + 7d + 30s: ");
 //Serial.print(future.year(), DEC);
 //Serial.print('/');
 //Serial.print(future.month(), DEC);
 //Serial.print('/');
 //Serial.print(future.day(), DEC);
 //Serial.print(' ');
 //Serial.print(future.hour(), DEC);
 //Serial.print(':');
 //Serial.print(future.minute(), DEC);
 //Serial.print(':');
 //Serial.print(future.second(), DEC);
 //Serial.println();
 //Serial.println();
 delay(3000);
 Serial.println("End line");
}

 void timeDisplay()
{
  DateTime t=rtc.now();
  lcd.clear();
  
  lcd.setCursor(3, 0);
  lcd.print(t.timestamp(DateTime::TIMESTAMP_TIME));

  lcd.setCursor(0, 1);
  // Display abbreviated Day-of-Week in the lower left corner
  lcd.print(t.dayOfTheWeek());

  lcd.setCursor(4,1);
  lcd.print(t.timestamp(DateTime::TIMESTAMP_DATE));
  
}
