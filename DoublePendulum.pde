Pendulum p;

void setup() {
    size(400, 400);
    p = new Pendulum(1000);
    p.setMass(20, 20);
    p.setLine(100, 100);
    p.setTheta(radians(145), radians(145));
    p.endSet();
    frameRate(30);
}

void draw() {
    background(255);
    translate(width / 2, height / 2);
    p.update();
    p.display();
}

void keyPressed() {
    if(key == 's'){
        save("frame/second.png");
    }
}