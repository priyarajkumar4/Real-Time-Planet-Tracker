

void getGPS() {


if (gps.location.isValid())
{ Lat = gps.location.lat();
Long = gps.location.lng();
}
else
{
Lat = 37.783393; //my location in case there is no gps readings
Long = -122.436845;

}

yy = gps.date.year();
//yy=2016;
//  Serial.println(yy);
mu = gps.date.month();
//mu=5;
// Serial.println(mu);
dd = gps.date.day();
// dd=4;
//  Serial.println(dd);
if (digitalRead(2) == 1)
{
if(gps.time.isValid())
{
hh = gps.time.hour();
//Serial.println(hh);
mm = gps.time.minute();
}
else
{
hh=3;
mm=0;
}
}
else if((digitalRead(2) == 0))
{
int ss=trajec();
hh = gps.time.hour();
mm = gps.time.minute();
}
}
