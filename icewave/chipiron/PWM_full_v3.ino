// =====================================================================
// PWM_full_v3.ino -- Chipiron Arduino Mega firmware (v3). /!\  LORA not tested.
//
// One firmware, two modes chosen at run time (no re-upload needed):
//   mc -> CONSTANT mode    accepts c<duty>; rest duty 10 %  (0 rpm with a UNIPOLAR ESCON config:
//                          10 % -> 0 rpm, 90 % -> +R rpm)
//   mo -> OSCILLATION mode accepts c<duty>, o<mean>,<amp>,<freq>; rest duty 50 %
//                          (0 rpm with the BIPOLAR ESCON config: 10 % -> -R, 50 % -> 0, 90 % -> +R rpm,
//                          e.g. parameters_Sep26.edc with R = 800)
// !! The ESCON configuration must match the mode, otherwise the rest duty does NOT give 0 rpm.
//
// START 
//     - START_ENABLED 0: at power-up / reset the ESCON enable line (DO1) is LOW (motor disabled, coasting);
//     - DEFAULT_MODE is MODE_OSCILLATION (bipolar, the current .edc);
//     - sending "mc" or "mo" selects the mode, writes its rest duty and THEN enables the motor (arming);
//     - "e" enables without changing mode (use only if you know the duty/mode are right).
//
// Commands (text, from USB serial = the phone, e.g. /usb-cmd/mo, and/or from LoRa):
//   mc | mo                 select mode (constant / oscillation), go to its rest duty, enable the motor
//   c<duty>                 constant duty in %, decimals accepted, e.g. c50  c60.5   (both modes)
//   o<mean>,<amp>,<freq>    sine: duty in %, freq in Hz, e.g. o50,20,0.5      (oscillation mode only)
//   x                       stop: back to the rest duty of the current mode (both modes)
//   e                       enable the ESCON (DO1 HIGH)
//   d                       disable the ESCON: duty back to rest, then DO1 LOW (motor coasts, NO braking:
//                           send x and wait for the ramp first if you want a controlled stop)
//   Examples:  mo   c60.0   o50.00,20.00,0.2000   x   d   e
//
// Telemetry (USB serial, 4 bytes every 250 ms):
//   b1 = ADC(A1) & 0xFF, b2 = ADC(A1) >> 4, b3 = ADC(A2) & 0xFF, b4 = ADC(A2) >> 4
//   decode: adc = (b2 << 4) | (b1 & 0x0F)     (see chipiron_decode_v3.py)
//   A1 = ESCON AO1 = actual CURRENT, A2 = ESCON AO2 = actual SPEED. The older comments said the opposite; the variable names telemSpeed / telemCurrent below
//   were NOT renamed and are therefore swapped: telemSpeed holds A1 (current), telemCurrent holds A2 (speed).
//
// Resolution: analogWrite is 8 bit -> duty steps of 100/255 = 0.39 % (~7.8 rpm with 20 rpm/%).
//
// Pins: PWM 13 -> ESCON digital input 1 (PWM set value), DO1 52 -> ESCON digital input 2 (enable),
//       DO2 51 -> direction (ignored by the ESCON as configured), A1 / A2 <- ESCON analog outputs.
// Used with: chipiron_session_v3.py, chipiron_run_v3.py, chipiron_oscillate_v3.py, mobile_v3.py (usb ...).
// =====================================================================

// ---------------- user switches ----------------
#define MODE_CONSTANT    0
#define MODE_OSCILLATION 1
#define DEFAULT_MODE MODE_OSCILLATION  // mode at power-up / reset: bipolar (current .edc); change with "mc" / "mo"
#define START_ENABLED 0     // 0 = ESCON disabled (DO1 LOW) until "mc" / "mo" / "e" is received; 1 = enabled at boot
#define USE_LORA 0          // 1 = also accept commands received by a Grove LoRa module on Serial1
#define LORA_ID_MOTOR  2    // this board's id in LoRa frames
#define LORA_ID_LAPTOP 1    // id of the laptop-side bridge (lora_tx_bridge.ino)
#define USE_ACK 1           // 1 = send an acknowledgement frame after each LoRa command
#define REST_DUTY_CONST 10.0   // constant mode: duty (%) giving 0 rpm for the unipolar mapping
#define REST_DUTY_OSC   50.0   // oscillation mode: duty (%) giving 0 rpm for bipolar config
// -----------------------------------------------

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#if USE_LORA
  #include <RH_RF95.h>                        // Grove LoRa library used in the BoxCo files
  RH_RF95<HardwareSerial> rf95(Serial1);      // Mega: Serial1 = pins 18 (TX1) / 19 (RX1)
  #define LORA_HEADER "HCHIP"                 // different from BoxCo's "HBXCO" on purpose
  #define LORA_ENDER  "ECHIP"
  #define LORA_FRAMES_BEFORE_REINIT 10        // BoxCo_Server re-inits the module every 10 frames
  uint8_t loraFrames = 0;
  bool loraOK = false;
  unsigned long lastLoraTry = 0;
#endif

#define PWM 13               // PWM output -> ESCON digital input 1 (PWM set value)
#define DO1 52               // ESCON enable (digital input 2)
#define DO2 51               // direction pin (ignored by ESCON as configured)

const int ai1 = A1;          // first analog telemetry channel (ESCON AO1 = actual current)
const int ai2 = A2;          // second analog telemetry channel (ESCON AO2 = actual speed)
// Wiring verified 2026-10-20: A1 <- AO1 (current), A2 <- AO2 (speed) (functions as set in the .edc; re-check if rewired).

const float DUTY_MIN = 10.0; // ESCON PWM input only accepts 10..90 %
const float DUTY_MAX = 90.0;

// ---- mode state ----
uint8_t mode = DEFAULT_MODE;
float restDuty = (DEFAULT_MODE == MODE_CONSTANT) ? REST_DUTY_CONST : REST_DUTY_OSC;
bool motorEnabled = (START_ENABLED != 0);   // "e" / "d" / "mc" / "mo" commands; DO1 follows this flag

// ---- wave state ----
float meanDuty = restDuty;   // %
float ampDuty  = 0.0;        // %
float freqHz   = 0.0;        // Hz
float phase    = 0.0;        // rad
unsigned long lastWave = 0;
const unsigned long WAVE_PERIOD_MS = 10;     // duty updated at 100 Hz

// ---- telemetry state machine (non-blocking timing) ----
const unsigned long TELEM_PERIOD_MS = 250;
unsigned long telemT0 = 0;
uint8_t telemStage = 0;
int telemSpeed = 0, telemCurrent = 0;

void writeDuty(float duty) {
  duty = constrain(duty, DUTY_MIN, DUTY_MAX);
  analogWrite(PWM, (int)(duty * 255.0 / 100.0 + 0.5));   // ~0.39 % resolution
}

void clampWave() {
  meanDuty = constrain(meanDuty, DUTY_MIN, DUTY_MAX);
  float room = min(meanDuty - DUTY_MIN, DUTY_MAX - meanDuty);
  ampDuty = constrain(ampDuty, 0.0, room);   // mean +/- amp must stay inside 10..90 %
  if (freqHz < 0) freqHz = 0;
}

// Select the mode and go back to its rest duty.
void setMode(uint8_t m) {
  mode = m;
  restDuty = (mode == MODE_CONSTANT) ? REST_DUTY_CONST : REST_DUTY_OSC;
  meanDuty = restDuty;
  ampDuty = 0.0;
  freqHz = 0.0;
  phase = 0.0;
  clampWave();
}

// Execute one command string such as "mo", "c50", "o50,20,0.5", "x", "e" or "d".
void applyCommand(char *s) {
  // Strip trailing newline or carriage return characters safely
  int len = strlen(s);
  while (len > 0 && (s[len - 1] == '\r' || s[len - 1] == '\n')) {
    s[len - 1] = '\0';
    len--;
  }

  // Tolerate URL-encoded commas ("%2C") if the phone server doesn't decode them
  char *p;
  while ((p = strstr(s, "%2C")) != NULL || (p = strstr(s, "%2c")) != NULL) {
    *p = ',';
    memmove(p + 1, p + 3, strlen(p + 3) + 1);
  }

  if (s[0] == 'm') {
    if (s[1] == 'c' || s[1] == 'o') {
      setMode(s[1] == 'c' ? MODE_CONSTANT : MODE_OSCILLATION);
      writeDuty(restDuty);               // rest duty on the PWM pin BEFORE the ESCON is enabled
      motorEnabled = true;               // arming: the mode is now known
    }
  }
  else if (s[0] == 'c') {
    meanDuty = atof(s + 1);
    ampDuty = 0.0; 
    freqHz = 0.0;
    clampWave();
  }
  else if (s[0] == 'o' && mode == MODE_OSCILLATION) {   // ignored in constant mode
    char *tok = strtok(s + 1, ",");
    if (tok) meanDuty = atof(tok);
    tok = strtok(NULL, ",");
    if (tok) ampDuty = atof(tok);
    tok = strtok(NULL, ",");
    if (tok) freqHz = atof(tok);
    phase = 0.0;                         // start at the mean value, going up
    clampWave();
  }
  else if (s[0] == 'x') {
    meanDuty = restDuty; 
    ampDuty = 0.0; 
    freqHz = 0.0;
  }
  else if (s[0] == 'e') {                // enable the motor
    motorEnabled = true;
  }
  else if (s[0] == 'd') {                // disable the motor, with the command back at rest
    meanDuty = restDuty;
    ampDuty = 0.0;
    freqHz = 0.0;
    motorEnabled = false;
  }
}

void readSerialCommands() {
  if (!Serial.available()) return;
  delay(5);                                   // let the rest of the message arrive
  char rec[32];
  int i = 0;
  while (Serial.available() && i < 31) rec[i++] = Serial.read();
  rec[i] = '\0';
  while (Serial.available()) Serial.read();  // drop excess buffer junk
  applyCommand(rec);
}

#if USE_LORA
void loraInit() {
  loraOK = rf95.init();
  if (loraOK) rf95.setFrequency(868.0);      
  lastLoraTry = millis();
}

void loraAck(const char *cmd) {
#if USE_ACK
  char m[10], a[10], f[10], frame[80];
  dtostrf(meanDuty, 1, 2, m);
  dtostrf(ampDuty, 1, 2, a);
  dtostrf(freqHz, 1, 4, f);
  snprintf(frame, sizeof(frame), "%s\t%d\t%d\t1\tOK %s m=%s a=%s f=%s\t%s",
           LORA_HEADER, LORA_ID_LAPTOP, LORA_ID_MOTOR, cmd, m, a, f, LORA_ENDER);
  rf95.send((uint8_t *)frame, strlen(frame));
  rf95.waitPacketSent();
#endif
}

void pollLora() {
  if (!loraOK) {                            
    if (millis() - lastLoraTry > 5000) loraInit();
    return;
  }
  if (!rf95.available()) return;
  uint8_t buf[RH_RF95_MAX_MESSAGE_LEN];
  uint8_t len = sizeof(buf);
  if (!rf95.recv(buf, &len)) return;
  if (len >= sizeof(buf)) len = sizeof(buf) - 1;
  buf[len] = '\0';

  char *tok = strtok((char *)buf, "\t");
  if (!tok || strcmp(tok, LORA_HEADER) != 0) return;
  tok = strtok(NULL, "\t");
  if (!tok || atoi(tok) != LORA_ID_MOTOR) return;       
  strtok(NULL, "\t");                                  
  strtok(NULL, "\t");                                  
  char *cmdTok = strtok(NULL, "\t");
  char *endTok = strtok(NULL, "\t");
  if (!cmdTok || !endTok || strcmp(endTok, LORA_ENDER) != 0) return;

  char cmd[32];
  strncpy(cmd, cmdTok, sizeof(cmd) - 1);
  cmd[sizeof(cmd) - 1] = '\0';
  applyCommand(cmd);
  loraAck(cmdTok);

  if (++loraFrames >= LORA_FRAMES_BEFORE_REINIT) {      
    loraFrames = 0;
    loraInit();
  }
}
#endif

void waveUpdate(unsigned long now) {
  if (now - lastWave < WAVE_PERIOD_MS) return;
  float dt = (now - lastWave) / 1000.0;
  lastWave = now;

  if (ampDuty > 0.0 && freqHz > 0.0) {
    phase += 2.0 * PI * freqHz * dt;
    while (phase >= 2.0 * PI) phase -= 2.0 * PI;
    writeDuty(meanDuty + ampDuty * sin(phase));
  } else {
    writeDuty(meanDuty);
  }
}

// Telemetry state machine (250 ms cycle)
void telemetry(unsigned long now) {
  unsigned long dt = now - telemT0;
  if (telemStage == 0) {
    telemSpeed = analogRead(ai1);                // NB: A1 = CURRENT (name kept, see header)
    Serial.write((uint8_t)(telemSpeed & 0xFF));
    telemStage = 1;
  } else if (telemStage == 1 && dt >= 10) {
    Serial.write((uint8_t)(telemSpeed >> 4));
    telemStage = 2;
  } else if (telemStage == 2 && dt >= 50) {
    telemCurrent = analogRead(ai2);              // NB: A2 = SPEED (name kept, see header)
    Serial.write((uint8_t)(telemCurrent & 0xFF));
    telemStage = 3;
  } else if (telemStage == 3 && dt >= 60) {
    Serial.write((uint8_t)(telemCurrent >> 4));
    telemStage = 4;
  } else if (telemStage == 4 && dt >= TELEM_PERIOD_MS) {
    telemT0 += TELEM_PERIOD_MS;
    telemStage = 0;
  }
}

void setup() {
  Serial.begin(115200);
  pinMode(PWM, OUTPUT);
  pinMode(DO1, OUTPUT);
  pinMode(DO2, OUTPUT);
  digitalWrite(DO1, motorEnabled ? HIGH : LOW);   // enable only if START_ENABLED
  digitalWrite(DO2, HIGH);
  writeDuty(restDuty);
#if USE_LORA
  loraInit();
#endif
}

void loop() {
  unsigned long now = millis();
  readSerialCommands();
#if USE_LORA
  pollLora();
#endif
  waveUpdate(now);
  telemetry(now);

  // Keep the ESCON enable line in the state requested by "e" / "d" (enabled by default)
  digitalWrite(DO1, motorEnabled ? HIGH : LOW);
  digitalWrite(DO2, HIGH);
}
