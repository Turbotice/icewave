// Pins
#define ENCA 13
#define ENCB 4

#define PWM 13
#define DO1 52  // to enable the motor, not use currently
#define DO2 51

int ai1 = A1; // pin speed ESCON 70/10 connected to analog pin 1
int ai2 = A2; // pin current ESCON 70/10 connected to analog pin 2



int speed = 0;  // variable to store the value read
int current = 0;  // variable to store the value read

int cmd = 10;

void setup() {
  Serial.begin(115200);

  pinMode(ENCA,INPUT);
  pinMode(ENCB,INPUT);
  pinMode(PWM,OUTPUT);
  pinMode(DO1,OUTPUT);
  pinMode(DO2,OUTPUT);

}

void loop() {
  char rec[10];// be careful with the size
  int i=0;
  int j=0;
  String s;

  speed = read_speed(ai1);
  Serial.write(speed);
  delay(10);
  Serial.write(speed/16);

  delay(40);
  current = read_current(ai2);
  Serial.write(current);
  delay(10);
  Serial.write(current/16);
  delay(140);

  while (Serial.available()) {        // If anything comes in Serial (USB),
    rec[i] = Serial.read();  // read it
    i++;
  }
  if (i>0){
    if (rec[0]=='c'){
    for (j=1;j<i;j++){
      //Serial.println(rec[j]);
      s = s+rec[j];
    }
    cmd = s.toInt();
    //Serial.println(cmd);
    }
      digitalWrite(DO1,LOW); 
    }

  delay(50);
  int pwr = int(cmd*256/100);
  //Serial.println(pwr);

  int dir = 1;
  setMotor(dir,pwr,PWM,DO1,DO2);

}

int read_speed(int ai1){
  speed = analogRead(ai1);  // read the input pin
  //Serial.println(speed); 
  return speed;      // debug value
}

int read_current(int ai2){
  current = analogRead(ai2);  // read the input pin
  //Serial.println(current);          // debug value
  return current;
}



void setMotor(int dir, int pwmVal, int pwm, int out1,int out2){
  if(dir == 1){ 
    // Turn one way
    digitalWrite(out1,HIGH); 
    digitalWrite(out2,HIGH);

     }
  else if(dir == -1){
    // Turn the other way
    digitalWrite(out1,HIGH);
    digitalWrite(out2,LOW);

  }
  else{
    // Or dont turn
    digitalWrite(out1,LOW);
    digitalWrite(out2,LOW);
  }
  analogWrite(pwm,pwmVal); // Motor speed
}
