String incomingByte ;    

const int ByPasson =  5;
const int ByPassoff =  10;
const int ledpin =  13;

void setup() {

  Serial.begin(9600);

  pinMode(ByPasson, OUTPUT);

  pinMode(ByPassoff, OUTPUT);

  pinMode(ledpin, OUTPUT);

}
void loop() {

  if (Serial.available() > 0) {

  incomingByte = Serial.readStringUntil('\n');

    if (incomingByte == "on") {

      digitalWrite(ByPasson, HIGH);

      delay(5000);

      digitalWrite(ByPassoff, HIGH);

      digitalWrite(ledpin, HIGH);

      Serial.write("ByPass on");

    }

    else if (incomingByte == "off") {

      digitalWrite(ByPassoff, LOW);

      delay(5000);

      digitalWrite(ByPasson, LOW);

      digitalWrite(ledpin, LOW);

      Serial.write("ByPass off");

    }


  }

}

