void MadgwickQuaternionUpdate(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
{
  float q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];   // short name local variable for readability
  float norm;
  float hx, hy, _2bx, _2bz;
  float s1, s2, s3, s4;
  float qDot1, qDot2, qDot3, qDot4;

  // Auxiliary variables to avoid repeated arithmetic
  float _2q1mx;
  float _2q1my;
  float _2q1mz;
  float _2q2mx;
  float _4bx;
  float _4bz;
  float _2q1 = 2.0f * q1;
  float _2q2 = 2.0f * q2;
  float _2q3 = 2.0f * q3;
  float _2q4 = 2.0f * q4;
  float _2q1q3 = 2.0f * q1 * q3;
  float _2q3q4 = 2.0f * q3 * q4;
  float q1q1 = q1 * q1;
  float q1q2 = q1 * q2;
  float q1q3 = q1 * q3;
  float q1q4 = q1 * q4;
  float q2q2 = q2 * q2;
  float q2q3 = q2 * q3;
  float q2q4 = q2 * q4;
  float q3q3 = q3 * q3;
  float q3q4 = q3 * q4;
  float q4q4 = q4 * q4;

  // Normalise accelerometer measurement
  norm = sqrt(ax * ax + ay * ay + az * az);
  if (norm == 0.0f) return; // handle NaN
  norm = 1.0f / norm;
  ax *= norm;
  ay *= norm;
  az *= norm;

  // Normalise magnetometer measurement
  norm = sqrt(mx * mx + my * my + mz * mz);
  if (norm == 0.0f) return; // handle NaN
  norm = 1.0f / norm;
  mx *= norm;
  my *= norm;
  mz *= norm;

  // Reference direction of Earth's magnetic field
  _2q1mx = 2.0f * q1 * mx;
  _2q1my = 2.0f * q1 * my;
  _2q1mz = 2.0f * q1 * mz;
  _2q2mx = 2.0f * q2 * mx;
  hx = mx * q1q1 - _2q1my * q4 + _2q1mz * q3 + mx * q2q2 + _2q2 * my * q3 + _2q2 * mz * q4 - mx * q3q3 - mx * q4q4;
  hy = _2q1mx * q4 + my * q1q1 - _2q1mz * q2 + _2q2mx * q3 - my * q2q2 + my * q3q3 + _2q3 * mz * q4 - my * q4q4;
  _2bx = sqrt(hx * hx + hy * hy);
  _2bz = -_2q1mx * q3 + _2q1my * q2 + mz * q1q1 + _2q2mx * q4 - mz * q2q2 + _2q3 * my * q4 - mz * q3q3 + mz * q4q4;
  _4bx = 2.0f * _2bx;
  _4bz = 2.0f * _2bz;

  // Gradient decent algorithm corrective step
  s1 = -_2q3 * (2.0f * q2q4 - _2q1q3 - ax) + _2q2 * (2.0f * q1q2 + _2q3q4 - ay) - _2bz * q3 * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (-_2bx * q4 + _2bz * q2) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + _2bx * q3 * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
  s2 = _2q4 * (2.0f * q2q4 - _2q1q3 - ax) + _2q1 * (2.0f * q1q2 + _2q3q4 - ay) - 4.0f * q2 * (1.0f - 2.0f * q2q2 - 2.0f * q3q3 - az) + _2bz * q4 * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (_2bx * q3 + _2bz * q1) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + (_2bx * q4 - _4bz * q2) * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
  s3 = -_2q1 * (2.0f * q2q4 - _2q1q3 - ax) + _2q4 * (2.0f * q1q2 + _2q3q4 - ay) - 4.0f * q3 * (1.0f - 2.0f * q2q2 - 2.0f * q3q3 - az) + (-_4bx * q3 - _2bz * q1) * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (_2bx * q2 + _2bz * q4) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + (_2bx * q1 - _4bz * q3) * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
  s4 = _2q2 * (2.0f * q2q4 - _2q1q3 - ax) + _2q3 * (2.0f * q1q2 + _2q3q4 - ay) + (-_4bx * q4 + _2bz * q2) * (_2bx * (0.5f - q3q3 - q4q4) + _2bz * (q2q4 - q1q3) - mx) + (-_2bx * q1 + _2bz * q3) * (_2bx * (q2q3 - q1q4) + _2bz * (q1q2 + q3q4) - my) + _2bx * q2 * (_2bx * (q1q3 + q2q4) + _2bz * (0.5f - q2q2 - q3q3) - mz);
  norm = sqrt(s1 * s1 + s2 * s2 + s3 * s3 + s4 * s4);    // normalise step magnitude
  norm = 1.0f / norm;
  s1 *= norm;
  s2 *= norm;
  s3 *= norm;
  s4 *= norm;

  // Compute rate of change of quaternion
  qDot1 = 0.5f * (-q2 * gx - q3 * gy - q4 * gz) - beta * s1;
  qDot2 = 0.5f * (q1 * gx + q3 * gz - q4 * gy) - beta * s2;
  qDot3 = 0.5f * (q1 * gy - q2 * gz + q4 * gx) - beta * s3;
  qDot4 = 0.5f * (q1 * gz + q2 * gy - q3 * gx) - beta * s4;

  // Integrate to yield quaternion
  q1 += qDot1 * deltat;
  q2 += qDot2 * deltat;
  q3 += qDot3 * deltat;
  q4 += qDot4 * deltat;
  norm = sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4);    // normalise quaternion
  norm = 1.0f / norm;
  q[0] = q1 * norm;
  q[1] = q2 * norm;
  q[2] = q3 * norm;
  q[3] = q4 * norm;

}



// Similar to Madgwick scheme but uses proportional and integral filtering on the error between estimated reference vectors and
// measured ones.
void MahonyQuaternionUpdate(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
{
  float q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];   // short name local variable for readability
  float norm;
  float hx, hy, bx, bz;
  float vx, vy, vz, wx, wy, wz;
  float ex, ey, ez;
  float pa, pb, pc;

  // Auxiliary variables to avoid repeated arithmetic
  float q1q1 = q1 * q1;
  float q1q2 = q1 * q2;
  float q1q3 = q1 * q3;
  float q1q4 = q1 * q4;
  float q2q2 = q2 * q2;
  float q2q3 = q2 * q3;
  float q2q4 = q2 * q4;
  float q3q3 = q3 * q3;
  float q3q4 = q3 * q4;
  float q4q4 = q4 * q4;

  // Normalise accelerometer measurement
  norm = sqrt(ax * ax + ay * ay + az * az);
  if (norm == 0.0f) return; // handle NaN
  norm = 1.0f / norm;        // use reciprocal for division
  ax *= norm;
  ay *= norm;
  az *= norm;

  // Normalise magnetometer measurement
  norm = sqrt(mx * mx + my * my + mz * mz);
  if (norm == 0.0f) return; // handle NaN
  norm = 1.0f / norm;        // use reciprocal for division
  mx *= norm;
  my *= norm;
  mz *= norm;

  // Reference direction of Earth's magnetic field
  hx = 2.0f * mx * (0.5f - q3q3 - q4q4) + 2.0f * my * (q2q3 - q1q4) + 2.0f * mz * (q2q4 + q1q3);
  hy = 2.0f * mx * (q2q3 + q1q4) + 2.0f * my * (0.5f - q2q2 - q4q4) + 2.0f * mz * (q3q4 - q1q2);
  bx = sqrt((hx * hx) + (hy * hy));
  bz = 2.0f * mx * (q2q4 - q1q3) + 2.0f * my * (q3q4 + q1q2) + 2.0f * mz * (0.5f - q2q2 - q3q3);

  // Estimated direction of gravity and magnetic field
  vx = 2.0f * (q2q4 - q1q3);
  vy = 2.0f * (q1q2 + q3q4);
  vz = q1q1 - q2q2 - q3q3 + q4q4;
  wx = 2.0f * bx * (0.5f - q3q3 - q4q4) + 2.0f * bz * (q2q4 - q1q3);
  wy = 2.0f * bx * (q2q3 - q1q4) + 2.0f * bz * (q1q2 + q3q4);
  wz = 2.0f * bx * (q1q3 + q2q4) + 2.0f * bz * (0.5f - q2q2 - q3q3);

  // Error is cross product between estimated direction and measured direction of gravity
  ex = (ay * vz - az * vy) + (my * wz - mz * wy);
  ey = (az * vx - ax * vz) + (mz * wx - mx * wz);
  ez = (ax * vy - ay * vx) + (mx * wy - my * wx);
  if (Ki > 0.0f)
  {
    eInt[0] += ex;      // accumulate integral error
    eInt[1] += ey;
    eInt[2] += ez;
  }
  else
  {
    eInt[0] = 0.0f;     // prevent integral wind up
    eInt[1] = 0.0f;
    eInt[2] = 0.0f;
  }

  // Apply feedback terms
  gx = gx + Kp * ex + Ki * eInt[0];
  gy = gy + Kp * ey + Ki * eInt[1];
  gz = gz + Kp * ez + Ki * eInt[2];

  // Integrate rate of change of quaternion
  pa = q2;
  pb = q3;
  pc = q4;
  q1 = q1 + (-q2 * gx - q3 * gy - q4 * gz) * (0.5f * deltat);
  q2 = pa + (q1 * gx + pb * gz - pc * gy) * (0.5f * deltat);
  q3 = pb + (q1 * gy - pa * gz + pc * gx) * (0.5f * deltat);
  q4 = pc + (q1 * gz + pa * gy - pb * gx) * (0.5f * deltat);

  // Normalise quaternion
  norm = sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4);
  norm = 1.0f / norm;
  q[0] = q1 * norm;
  q[1] = q2 * norm;
  q[2] = q3 * norm;
  q[3] = q4 * norm;

}
static void smartDelay(unsigned long ms)
{
  //Serial.println("C");
  unsigned long start = millis();
  do
  {
    while (Serial2.available())
      gps.encode(Serial2.read());
  } while (millis() - start < ms);
  //Serial.println("D");
}


// Similar to Madgwick scheme but uses proportional and integral filtering on the error between estimated reference vectors and
// measured ones.
int planetInput(int potval)
  
{

//int potval = analogRead(3); 


  if ((potval >= 0) && (potval <= 90))
  {
    Serial.println("Mercury");
    return 1;
  }
  else if ((potval >= 91) && (potval <= 191))
  {
    Serial.println("Venus");
    return 2;
  }
  else if ((potval >= 192) && (potval <= 292))
  {
    Serial.println("Mars");
    return 4;
  }
  else if ((potval >= 293) && (potval <= 393))
  {
    Serial.println("Jupiter");
    return 5;
  }
  else if ((potval >= 394) && (potval <= 494))
  {
    Serial.println("Saturn");
    return 6;
  }
  else if ((potval >= 495) && (potval <= 595))
  {
    Serial.println("Uranus");
    return 7;
  }
  else if ((potval >= 596) && (potval <= 696))
  {
    Serial.println("Neptune");
    return 8;
  }
  else if ((potval >= 697) && (potval <= 1023))
  {
    Serial.println("Pluto");
    return 9;
  }

//  if ((potval <= 700))
//  {
//    Serial.println("Mercury");
//    return 1;
//  }
//  else if ((potval >= 701))
//  {
//    Serial.println("Neptune");
//    return 8;
//  }
//  else if ((potval >= 801))
//  {
//    Serial.println("Mars");
//    return 4;
//  }
//  else if ((potval >= 601) && (potval <= 800))
//  {
//    Serial.println("Jupiter");
//    return 5;
//  }
//  else if ((potval >= 801) && (potval <= 1023))
//  {
//    Serial.println("Saturn");
//    return 6;
//  }
//  else if ((potval >= 495) && (potval <= 595))
//  {
//    Serial.println("Uranus");
//    return 7;
//  }
//  else if ((potval >= 596) && (potval <= 696))
//  {
//    Serial.println("Neptune");
//    return 8;
//  }
//  else if ((potval >= 697) && (potval <= 1023))
//  {
//    Serial.println("Pluto");
//    return 9;
//  }
}


int trajec()
{
  hh=0;
  mm=0;
 while(hh<24&& (digitalRead(46) == 0))
 {
   while(mm<60&& (digitalRead(46) == 0))
   {
     
mainCalculations();


nyaw = 360 - yaw;

  Azim = Azimuth - nyaw;
  Azim -= 90;
  while (Azim < 0)
    Azim = 360.0 - abs(Azim);

  Azi = map(Azim, 0, 360, 5, 29);
  Az = (int)Azi;
  Elev = map(Elevation, -90, 90, 2, 178);
  El = (int)Elev;
  myservoAz.write(Az);
  myservoEl.write(El);
  mm++;
  Serial.print(hh); Serial.print(":"); Serial.println(mm);
   }
   mm=0;
   hh++;
   while(hh==19)
   hh=20;
   Serial.println("good");
 }
 return 5;
}

double ipart(double xx)
{
  //Serial.println("IPART");
  double sgn;
  if (xx < 0)
  {
    sgn = -1.0;
  }

  else if (xx == 0)
  {
    sgn = 0.0;
  }

  else if (xx > 0)
  {
    sgn = 1.0;
  }
  double ret = sgn * ((int)fabs(xx));

  return ret;
}


double FNdegmin(double xx)
{
  //Serial.println("DEGMIN");
  double a = ipart(xx) ;
  double b = xx - a ;
  double e = ipart(60 * b) ;
  //   deal with carry on minutes
  if ( e >= 60 )
  {
    e = 0 ;
    a = a + 1 ;
  }
  return (a + (e / 100) );
}

double dayno(int dx, int mx, int yx, double fx)
{
//Serial.println("DAY NO");

  //dno=(367 * yx) -  (int)(7*(yx + (int)((mx + 9) / 12)) / 4) + (int)(275 * mx / 9) + dx - 730531.5 + fx;
  dno = 987 + dx + (fx / 24);
  //Serial.print("\ndays:");
  //Serial.println(dno);
  return dno;
}
double frange(double x)
{
  //Serial.println("FRANGE");
  x = x / (2 * pi);
  x = (2 * pi) * (x - ipart(x));
  if (x < 0)
    x = x + (2 * pi);
  return x;
}
double fkep( double m, double ecc)
{
 
    //  Serial.println(m);  Serial.println(ecc);
  //Serial.println("FKEPA");
 // Serial.println(m);  Serial.println(ecc);
   double e = ecc;

//  do
//  { dyo = e - (ecc * sin(e)) - m;
//    e = e - (dyo / (1 - (ecc * cos(e))));
//    delay(1);
////Serial.print("fabs");Serial.print((fabs(dyo)*pow(10, 6)+1),8);Serial.print("              ");Serial.print("pow");Serial.println((pow(10, -6)+1),8);
//  } while (fabs(dyo)>= pow(10, -12));
//      //Serial.println("fkepB");
//   double v = 2 * atan(sqrt((1 + ecc) / (1 - ecc)) * tan(e / 2));
   
 double v = m + (2 * e - 0.25 *pow(e,3) + 5/96 * pow(e,5)) * sin(m) + (1.25 * pow(e,2) - 11/24 * pow(e,4)) * sin(2*m) + (13/12 * pow(e,3) - 43/64 * pow(e,5)) * sin(3*m) + 103/96 * pow(e,4) * sin(4*m) + 1097/960 * pow(e,5) * sin(5*m);

  if (v < 0)
    v = v + (2 * pi);
  return v;
}
double fnatan(double x, double y)
{
  //Serial.println("ATAN");
  double a = atan(y / x);
  if (x < 0)
    a = a + pi;
  if ((y < 0) && (x > 0))
    a = a + (2 * pi);
  return a;
}
void AltAzCalculate(double RA, double Dec, double Lat, double Long, double hrs, double minut, double dy)
{
  //Serial.println("G");
  // Day offset and Local Siderial Time
  dy = dy + 4975.5;

  double LST = (100.46 + 0.985647 * dy + Long + 15 * (hrs + minut / 60) + 360) - (((int)((100.46 + 0.985647 * dy + Long + 15 * (hrs + minut / 60) + 360) / 360)) * 360);

  // Hour Angle
  double HA = (LST - RA + 360) - ((int)((LST - RA + 360) / 360)) * 360 ;

  // HA, DEC, Lat to Alt, AZ
  double x = cos(HA * (pi / 180)) * cos(Dec * (pi / 180));
  double y = sin(HA * (pi / 180)) * cos(Dec * (pi / 180));
  double z = sin(Dec * (pi / 180));

  double xhor = x * cos((90 - Lat) * (pi / 180)) - z * sin((90 - Lat) * (pi / 180));
  double yhor = y;
  double zhor = x * sin((90 - Lat) * (pi / 180)) + z * cos((90 - Lat) * (pi / 180));

  Azimuth = atan2(yhor, xhor) * (180 / pi) + 180;
  Elevation = asin(zhor) * (180 / pi);
}
void earth()
{
//Serial.println("B");
  M[3] = ((n[3] * rads) * d) + (L[3] - p[3]) * rads;
  M[3] = frange(M[3]);
  v[3] = fkep(M[3], e[3]);
  r[3] = a[3] * ((1 - (pow(e[3], 2))) / (1 + (e[3] * cos(v[3]))));
  x[3] = r[3] * cos(v[3] + p[3] * rads);
  y[3] = r[3] * sin(v[3] + p[3] * rads);
  z[3] = 0;
}
void mainCalculations()
{

  dfrac = hh + (mm / 60);
 d = dayno(dd, mu, yy, dfrac);
 //Serial.println("E");
 earth();
 //Serial.println("F");
//Serial.println("E");
 int j;
 for (j = 0; j <=9; j++)
  {
  if (j == 3)
      continue;
      if(j==9);
  // Serial.println("A");
    M[j] = ((n[j] * rads) * d) + (L[j] - p[j]) * rads;
      if(j==9);
   //   Serial.println("B");

   M[j] = frange(M[j]);
     if(j==9);
   // Serial.println("C");
    v[j] = fkep(M[j], e[j]);
      if(j==9);
     // Serial.println("D");

   r[j] = a[j] * ((1 - pow(e[j], 2)) / (1 + (e[j] * cos(v[j]))));
   x[j] = r[j] * (cos(o[j] * rads) * cos(v[j] + p[j] * rads - o[j] * rads) - sin(o[j] * rads) * sin(v[j] + p[j] * rads - o[j] * rads) * cos(i[j] * rads));
   y[j] = r[j] * (sin(o[j] * rads) * cos(v[j] + p[j] * rads - o[j] * rads) + cos(o[j] * rads) * sin(v[j] + p[j] * rads - o[j] * rads) * cos(i[j] * rads));
   z[j] = r[j] * (sin(v[j] + p[j] * rads - o[j] * rads)) * sin(i[j] * rads);
   Xi[j] = x[j] - x[3];
   Yi[j] = y[j] - y[3];
   Zi[j] = z[j];

   Xq[j] = Xi[j];
    Yq[j] = (Yi[j] * cos(ec)) - (Zi[j] * sin(ec));
    Zq[j] = (Yi[j] * sin(ec)) + (Zi[j] * cos(ec));
    ra[j] = fnatan(Xq[j], Yq[j]);
    dec[j] = atan(Zq[j] / sqrt(pow(Xq[j], 2.0) + pow(Yq[j], 2.0)));
  // Serial.println(j);
  }

//Serial.println("H");
  double alpha = FNdegmin((ra[pno] * degs) / 15);
  double delta = FNdegmin(dec[pno] * degs);

//Serial.println("G");


  AltAzCalculate((alpha * 15), delta, Lat, Long, hh, mm, d);
}
