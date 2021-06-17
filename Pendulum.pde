class Pendulum {
    // 試行回数
    int n = 0;
    int nowN = 0;

    // deltaTime
    float dTime = 0.01;

    // 質点の質量
    float mass1 = 0;
    float mass2 = 0;

    // 質点の座標
    float x1 = 0;
    float y1 = 0;
    float x2 = 0;
    float y2 = 0;

    // 線の長さ
    float line1 = 0;
    float line2 = 0;

    // 角度
    float theta1 = 0;
    float theta2 = 0;

    // 角速度
    float omega1 = 0;
    float omega2 = 0;

    // 重力加速度
    float gravity = 9.8;

    // μ
    float mu;

    // dθ
    float dTheta = theta1 - theta2;
    float cosTheta = cos(dTheta);
    float sinTheta = sin(dTheta);

    Pendulum(int _n){
        this.n = _n;
    }

    void endSet(){
        mu = 1 + (mass1 / mass2);
        dTheta = theta1 - theta2;
        cosTheta = cos(dTheta);
        sinTheta = sin(dTheta);
    }

    void update(){
        rungeKutta();
    }

    void display(){
        x1 = int(line1 * sin(theta1));
        x2 = int(x1 + line2 * sin(theta2));
        y1 = int(line1 * cos(theta1));
        y2 = int(y1 + line2 * cos(theta2));
        line(0, 0, x1, y1);
        line(x1, y1, x2, y2);
        fill(255, 0, 0);
        ellipse(x1, y1, 10, 10);
        fill(0, 0, 255);
        ellipse(x2, y2, 10, 10);

        this.nowN++;
        if(this.nowN < this.n) rec();
    }

    void rec(){
        println("N=" + this.nowN + " x1: " + x1 + " y1: " + y1 + " x2: " + x2 + " y2: " + y2);
    }

    // xn, ynを入れる配列
    private float[][] j = new float[2][4];
    // kn(傾き)を入れる配列
    private float[][] k = new float[2][4];

    // ルンゲ＝クッタ法による角速度の計算
    void rungeKutta(){
        j[0][0] = this.omega1;
        j[1][0] = this.omega2;
        k[0][0] = dTime * (gravity * (sin(this.theta2) * cosTheta - this.mu * sin(this.theta1))
          - (this.line1 * pow(this.omega1, 2) * cosTheta + this.line2 * pow(this.omega2, 2)) * sinTheta) / (this.line1 * (this.mu - pow(cosTheta, 2)));
        k[1][0] = dTime * (gravity * this.mu * (sin(this.theta1) * cosTheta - sin(this.theta2))
          - (this.mu * this.line1 * pow(this.omega1, 2) + this.line2 * pow(this.omega2, 2) * cosTheta) * sinTheta) / (this.line2 * (this.mu - pow(cosTheta, 2)));

        dTheta = (theta1 + j[0][0] * 0.5) - (theta2 + j[1][0] * 0.5);
        cosTheta = cos(dTheta);
        sinTheta = sin(dTheta);

        j[0][1] = this.omega1 + k[0][0] * 0.5;
        j[1][0] = this.omega2 + k[1][0] * 0.5;
        k[0][1] = rk1(1);
        k[1][1] = rk2(1);

        dTheta = (theta1 + j[0][1] * 0.5) - (theta2 + j[1][1] * 0.5);
        cosTheta = cos(dTheta);
        sinTheta = sin(dTheta);

        j[0][2] = this.omega1 + k[0][1] * 0.5;
        j[1][2] = this.omega2 + k[1][1] * 0.5;
        k[0][2] = rk1(2);
        k[1][2] = rk2(2);

        dTheta = (theta1 + j[0][1] * 0.5) - (theta2 + j[1][1] * 0.5);
        cosTheta = cos(dTheta);
        sinTheta = sin(dTheta);

        j[0][3] = this.omega1 + k[0][2] * 0.5;
        j[1][3] = this.omega2 + k[1][2] * 0.5;
        k[0][3] = rk1(3);
        k[1][3] = rk2(3);

        this.theta1 += (j[0][1] + (2 * j[0][1]) + (2 * j[0][2]) + j[0][3]) / 6;
        this.theta2 += (j[1][1] + (2 * j[1][1]) + (2 * j[1][2]) + j[1][3]) / 6;
        this.omega1 += (k[0][1] + (2 * k[0][1]) + (2 * k[0][2]) + k[0][3]) / 6;
        this.omega2 += (k[1][1] + (2 * k[1][1]) + (2 * k[1][2]) + k[1][3]) / 6;
    }

    float rk1(int i){
        return dTime * (gravity * (sin(this.theta2 + j[1][i-1] * 0.5) * cosTheta - this.mu * sin(this.theta1 + j[0][i-1] * 0.5)) 
            - (this.line1 * pow(this.omega1 + k[0][i-1] * 0.5, 2) * cosTheta + this.line2 * pow(this.omega2 + k[1][i-1] * 0.5, 2)) * sinTheta) / (this.line1 * (this.mu - pow(cosTheta, 2)));
    }

    float rk2(int i){
        return dTime * (gravity * this.mu * (sin(this.theta1 + j[1][i-1] * 0.5) * cosTheta - sin(this.theta2 + j[1][i-1] * 0.5))
            - (this.mu * this.line1 * pow(this.omega1 + k[0][i-1] * 0.5, 2) + this.line2 * pow(this.omega2 + k[1][i-1] * 0.5, 2) * cosTheta) * sinTheta) / (this.line2 * (this.mu - pow(cosTheta, 2)));
    }

    public void setMass1(float _m1){
        this.mass1 = _m1;
    }

    public void setMass2(float _m2){
        this.mass2 = _m2;
    }

    public void setMass(float _m1, float _m2){
        this.mass1 = _m1;
        this.mass2 = _m2;
    }

    public void setLine1(float _l1){
        this.line1 = _l1;
    }

    public void setLine2(float _l2){
        this.line2 = _l2;
    }

    public void setLine(float _l1, float _l2){
        this.line1 = _l1;
        this.line2 = _l2;
    }

    public void setTheta1(float _theta1){
        this.theta1 = _theta1;
    }

    public void setTheta2(float _theta2){
        this.theta2 = _theta2;
    }

    public void setTheta(float _theta1, float _theta2){
        this.theta1 = _theta1;
        this.theta2 = _theta2;
    }
}
