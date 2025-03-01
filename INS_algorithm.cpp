#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;

bool alignment_flag = 0;
bool flag = 0;

const float g = 9.81;
//shirota & dolgota / f % a
float latitude_0;
float longitude_0;

float R_Earth = 6.371 * pow(10, 6);
float R_latitude;  
float  R_longitude;   

float C_0;

float wX, wY, wZ;
float aX, aY, aZ;

float PITCH_0 = 0;
float ROLL_0 = 0;
float YAW_0 = 0;

float PITCH = 0;
float ROLL = 0;
float YAW = 0;

float wE, wN, wUp;

float current_aX = 0;
float previous_aX = 0;

float current_aY = 0; 
float previous_aY = 0;

float current_aZ = 0; 
float previous_aZ = 0;

float velocityX = 0;
float prev_velocityX = 0;

float velocityY = 0;
float prev_velocityY = 0;

float velocityZ = 0; 
float prev_velocityZ = 0;

float X, Y, Z;

float matrix_LL[3][3]; //матрица для выставки на первом такте 
float Acc_matrix_BL[3][1]; //матрица с показаниями акселерометров
float Gyro_matrix_BL[3][1]; //матрица с покзааниями ДУСов
float matrix_W_B[3][3]; //3x3 матрица показаний ДУСов
float matrix_W_LL[3][3]; //матрица посчитанных угловых скоростей из показаний акслерометра

float matrix(int i, int j){
    matrix_LL[0][0] = cos(ROLL) * cos(YAW) + sin(PITCH) * sin(ROLL) * sin(YAW);
    matrix_LL[0][1] = -cos(ROLL) * sin(YAW) + sin(PITCH) * sin(ROLL) * cos(YAW);
    matrix_LL[0][2] = cos(PITCH) * sin(ROLL);
    matrix_LL[1][0] = cos(PITCH) * sin(YAW);
    matrix_LL[1][1] = cos(PITCH) * cos(YAW);
    matrix_LL[1][2] = sin(PITCH);
    matrix_LL[2][0] = sin(ROLL) * cos(YAW) - sin(PITCH) * cos(ROLL) * sin(YAW);
    matrix_LL[2][1] = -sin(ROLL) * sin(YAW) - sin(PITCH) * cos(ROLL) * sin(YAW);
    matrix_LL[2][2] = cos(PITCH) * cos(ROLL);
    return matrix_LL[i][j];
}

float matrix_W_DUS(int i, int j){
    matrix_W_B[0][0] = 0;
    matrix_W_B[0][1] = -Gyro_matrix_BL[2][0];
    matrix_W_B[0][2] = Gyro_matrix_BL[1][0];
    matrix_W_B[1][0] = Gyro_matrix_BL[2][0];
    matrix_W_B[1][1] = 0;
    matrix_W_B[1][2] = -Gyro_matrix_BL[0][0];
    matrix_W_B[2][0] = -Gyro_matrix_BL[1][0];
    matrix_W_B[2][1] = Gyro_matrix_BL[0][0];
    matrix_W_B[2][2] = 0;
    return matrix_W_B[i][j];
}

float matrix_W_ENUp(int i, int j){
    matrix_W_LL[0][0] = 0;
    matrix_W_LL[0][1] = -wUp;
    matrix_W_LL[0][2] = wN;
    matrix_W_LL[1][0] = wUp;
    matrix_W_LL[1][1] = 0;
    matrix_W_LL[1][2] = -wE;
    matrix_W_LL[2][0] = -wN;
    matrix_W_LL[2][1] = wE;
    matrix_W_LL[2][2] = 0;
    return matrix_W_LL[i][j];
}

float bodyToLocal(float aX, float aY, float aZ, int i, int j){
    Acc_matrix_BL[0][0] = aX;
    Acc_matrix_BL[1][0] = aY;
    Acc_matrix_BL[2][0] = aZ;
    float Acc_matrix_ENUp[3][1] = {0, 0, 0}; //матрица ускорений в системе ENUp
    for(int i = 0, j = 0, k = 0; j <= 2; j++){
        while(k <= 2){
            Acc_matrix_ENUp[j][i] += ((matrix(j, k) * Acc_matrix_BL[k][i]));
            k++;
        }
        k = 0;
    }
    return Acc_matrix_ENUp[j][i];
}

float getX(float accel){
    previous_aX = current_aX;
    current_aX = accel;
    prev_velocityX = velocityX;
    velocityX += (((current_aX + previous_aX) / 2 ) * 0.005);
    //checkout logic of math operations below - divide and integrate
    X = latitude_0;
    X += ((((velocityX + prev_velocityX)  / 2 ) * 0.005) / R_latitude);
    return X;
}
  
  float getY(float accel){
    previous_aY = current_aY;
    current_aY = accel;
    prev_velocityY = velocityY;
    velocityY += (((current_aY + previous_aY) / 2 ) * 0.005);
    //checkout logic of math operations below - divide and integrate
    Y = longitude_0;
    Y += ((((velocityY + prev_velocityY) / 2 ) * 0.005) / R_longitude);
    return Y;
}

float Poisson(int i, int j){
    float NEW_LL_MATRIX[3][3] = {0}; // новая матрица, которая используется для перевода в систему ENUp на последующих тактах работы алгоритма
    for(int i = 0, j = 0, k = 0; i <= 2; i++){
      for(; j <= 2; j++){
        for(; k <= 2; k++){
          NEW_LL_MATRIX[i][j] += ((matrix(i, k) * matrix_W_DUS(k, j)) - matrix(i, k) * matrix_W_ENUp(k ,j));
        }
        k = 0;
      }
      j = 0;
    }
    return NEW_LL_MATRIX[i][j];
}

int main()
{
    float latitude = 0;
    //setlocale(LC_ALL, "ru");
    // Путь к документу
    //string path = "test_data.txt";
    string path = "20170922_Cessna172_200Hz_Ref.txt";
    ifstream Fin;
    Fin.open(path);

    if (!Fin.is_open()) {
        cerr << "Ошибка открытия файла!" << endl;
        return 1;
    }

    string str; // Переменная для записи строк
    int numLines = 0;

    while (getline(Fin, str)) {
        ++numLines;

        // Разбиение на строки при помощи std::istringstream
        istringstream iss(str);
        string token;
        int tokenCount = 0;

        while (iss >> token) {
            ++tokenCount;
            if(tokenCount == 2){
                //wX
                Gyro_matrix_BL[0][0] = stof(token);
            }

            if(tokenCount == 3){
                //wY
                Gyro_matrix_BL[1][0] = stof(token);
            }

            if(tokenCount == 4){
                //wZ
                Gyro_matrix_BL[2][0] = stof(token);
            }

            if(tokenCount == 5){
                //aX
                Acc_matrix_BL[0][0] = stof(token);             
            }

            if(tokenCount == 6){
                //aY
                Acc_matrix_BL[1][0] = stof(token);
            }

            if(tokenCount == 7){
                //aZ
                Acc_matrix_BL[2][0] = stof(token);   
            }
            if(tokenCount == 13){
                //shirota
                latitude_0 = stof(token);
            }
        }

        //cout << numLines << "    " << "omegaX = " << Gyro_matrix_BL[0][0] << "   " << "omegaY = " << Gyro_matrix_BL[1][0] << "   " << "omegaZ = " << Gyro_matrix_BL[2][0] << endl;
        //alignment
        if(Gyro_matrix_BL[0][0] == 0){
            alignment_flag = 1;
            PITCH_0 = -Acc_matrix_BL[0][0] / g;
            ROLL_0 = -Acc_matrix_BL[1][0] / g;
            YAW_0 = -atan(Gyro_matrix_BL[0][0] / Gyro_matrix_BL[1][0]);
            cout << PITCH_0 << "    " << ROLL_0 << "    " << YAW_0 << endl;
        }

        if(Gyro_matrix_BL[0][0] != 0){
            if(flag == 0){
                flag = 1;
                PITCH = PITCH_0;
                ROLL = ROLL_0;
                YAW = YAW_0;
                latitude = latitude_0;
                cout << PITCH << "    " << ROLL << "    " << YAW << endl;
            }
            else{

                //cout << wE << "    " << wN << "    " << wUp << endl;

                bodyToLocal(Acc_matrix_BL[0][0], Acc_matrix_BL[1][0], Acc_matrix_BL[2][0], 0, 0);
                

                X = getX(bodyToLocal(Acc_matrix_BL[0][0], Acc_matrix_BL[1][0], Acc_matrix_BL[2][0], 0, 0));
                Y = getY(bodyToLocal(Acc_matrix_BL[0][0], Acc_matrix_BL[1][0], Acc_matrix_BL[2][0], 1, 0));
                //Z = getZ(Acc_matrix_ENUp[2][0]);

                R_latitude = (R_Earth * (1 - pow(M_E, 2))) / pow(1 - pow(M_E, 2) * pow(sin(latitude_0), 2), 3 / 2);
                R_longitude = (R_Earth * (1 - pow(M_E, 2))) / pow(1 - pow(M_E, 2) * pow(sin(latitude_0), 2), 1 / 2);


                //calculate angular velocity
                wE = -velocityY / R_latitude;
                wN = velocityX / R_longitude;
                wUp = velocityX / R_longitude * tan(latitude_0);

                Poisson(0 ,0);

                C_0 = sqrt(((Poisson(2, 0) * Poisson(2, 0)) + (Poisson(2, 2) * Poisson(2, 2))));

                PITCH = atan(Poisson(2, 1) / C_0);//(NEW_LL_MATRIX[2][1] / C_0);
                ROLL = atan2(Poisson(2, 0), Poisson(2, 2));//(NEW_LL_MATRIX[2][0] / NEW_LL_MATRIX[2][2]);
                YAW = atan2(Poisson(0, 1), Poisson(1, 1)); //fmod(atan2(Poisson(0, 1), Poisson(1, 1)), 360);//((atan(fmod((Poisson(0, 1) / Poisson(1, 1)), 360.0));//(NEW_LL_MATRIX[0][1] / NEW_LL_MATRIX[1][1]);

                cout << PITCH << "    " << ROLL << "    " << YAW << endl;
            }
        }
    }

    Fin.close();

    cout << "Общее количество строк: " << numLines << endl;

    return 0;
}