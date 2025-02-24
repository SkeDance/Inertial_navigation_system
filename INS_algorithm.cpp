#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>

using namespace std;

const float g = 9.81;
//shirota & dolgota / f % a
float latitude_0 = 55.44;
float longitude_0 = 37.36;

float R_Earth = 6.371 * pow(10, 6);
float R_latitude = (R_Earth * (1 - pow(M_E, 2))) / pow(1 - pow(M_E, 2) * pow(sin(latitude_0), 2), 3 / 2);
float  R_longitude = (R_Earth * (1 - pow(M_E, 2))) / pow(1 - pow(M_E, 2) * pow(sin(longitude_0), 2), 1 / 2);

float C_0;

float wX, wY, wZ;
float aX, aY, aZ;
float PITCH_0, ROLL_0, YAW_0;
float PITCH, ROLL, YAW;

float wE, wN, wUp;

float current_aX, previous_aX;
float current_aY, previous_aY;
float current_aZ, previous_aZ;

float velocityX, prev_velocityX;
float velocityY, prev_velocityY;
float velocityZ, prev_velocityZ;

float X, Y, Z;

float matrix_LL[3][3]; //матрица для выставки на первом такте 
float Acc_matrix_BL[3][1]; //матрица с показаниями акселерометров
float Gyro_matrix_BL[3][1]; //матрица с покзааниями ДУСов
float Acc_matrix_ENUp[3][1] = {0, 0, 0}; //матрица ускорений в системе ENUp
float NEW_LL_MATRIX[3][3]; // новая матрица, которая используется для перевода в систему ENUp на последующих тактах работы алгоритма
float matrix_W_B[3][3]; //3x3 матрица показаний ДУСов
float matrix_W_LL[3][3]; //матрица посчитанных угловых скоростей из показаний акслерометра

float matrix(int i, int j){
    matrix_LL[0][0] = cos(ROLL_0) * cos(YAW_0) + sin(PITCH_0) * sin(ROLL_0) * sin(YAW_0);
    matrix_LL[0][1] = -cos(ROLL_0) * sin(YAW_0) + sin(PITCH_0) * sin(ROLL_0) * cos(YAW_0);
    matrix_LL[0][2] = cos(PITCH_0) * sin(ROLL_0);
    matrix_LL[1][0] = cos(PITCH_0) * sin(YAW_0);
    matrix_LL[1][1] = cos(PITCH_0) * cos(YAW_0);
    matrix_LL[1][2] = sin(PITCH_0);
    matrix_LL[2][0] = sin(ROLL_0) * cos(YAW_0) - sin(PITCH_0) * cos(ROLL_0) * sin(YAW_0);
    matrix_LL[2][1] = -sin(ROLL_0) * sin(YAW_0) - sin(PITCH_0) * cos(ROLL_0) * sin(YAW_0);
    matrix_LL[2][2] = cos(PITCH_0) * cos(ROLL_0);
    return matrix_LL[i][j];
}

float matrix_W_DUS(int i, int j){
    matrix_W_B[0][0] = 0;
    matrix_W_B[0][1] = -wZ;
    matrix_W_B[0][2] = wY;
    matrix_W_B[1][0] = wZ;
    matrix_W_B[1][1] = 0;
    matrix_W_B[1][2] = -wX;
    matrix_W_B[2][0] = -wY;
    matrix_W_B[2][1] = wX;
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

void bodyToLocal(float aX, float aY, float aZ){
    Acc_matrix_BL[0][0] = aX;
    Acc_matrix_BL[1][0] = aY;
    Acc_matrix_BL[2][0] = aZ;
    for(int i = 0, j = 0, k = 0; j <= 2; j++){
        while(k <= 2){
            Acc_matrix_ENUp[j][i] += ((matrix(j, k) * Acc_matrix_BL[k][i]));
            k++;
        }
        k = 0;
    }
}

float getX(float accel){
    previous_aX = current_aX;
    current_aX = accel;
    prev_velocityX = velocityX;
    velocityX += (((current_aX + previous_aX) / 2 ) * 0.01);
    //checkout logic of math operations below - divide and integrate
    X = latitude_0;
    X += ((((velocityX + prev_velocityX)  / 2 ) * 0.01) / R_latitude);
    return X;
}
  
  float getY(float accel){
    previous_aY = current_aY;
    current_aY = accel;
    prev_velocityY = velocityY;
    velocityY += (((current_aY + previous_aY) / 2 ) * 0.01);
    //checkout logic of math operations below - divide and integrate
    Y = longitude_0;
    Y += ((((velocityY + prev_velocityY) / 2 ) * 0.01) / R_longitude);
    return Y;
}

void Poisson(){
    for(int i = 0, j = 0, k = 0; i <= 2; i++){
      for(; j <= 2; j++){
        for(; k <= 2; k++){
          NEW_LL_MATRIX[i][j] += ((matrix(i, k) * matrix_W_DUS(k, j)) - matrix(i, k) * matrix_W_ENUp(k ,j));
        }
        k = 0;
      }
      j = 0;
    }
}

int main()
{
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
        }

        //cout << numLines << "    " << "omegaX = " << Gyro_matrix_BL[0][0] << "   " << "omegaY = " << Gyro_matrix_BL[1][0] << "   " << "omegaZ = " << Gyro_matrix_BL[2][0] << endl;

        if(numLines == 1){
            PITCH_0 = -Acc_matrix_BL[0][0] / g;
            ROLL_0 = -Acc_matrix_BL[1][0] / g;
            YAW_0 = -atan(Gyro_matrix_BL[0][0] / Gyro_matrix_BL[1][0]);
            cout << PITCH_0 << "    " << ROLL_0 << "    " << YAW_0 << endl;

            bodyToLocal(Acc_matrix_BL[0][0], Acc_matrix_BL[1][0], Acc_matrix_BL[2][0]);
            

            X = getX(Acc_matrix_ENUp[0][0]);
            Y = getY(Acc_matrix_ENUp[1][0]);
            //Z = getZ(Acc_matrix_ENUp[2][0]);

            //calculate angular velocity
            wE = -velocityY / R_latitude;
            wN = velocityX / R_longitude;
            wUp = velocityX / R_longitude * tan(latitude_0);

            Poisson();

            C_0 = sqrt(((NEW_LL_MATRIX[2][0] * NEW_LL_MATRIX[2][0]) + (NEW_LL_MATRIX[2][2] * NEW_LL_MATRIX[2][2])));

            PITCH = atan(fmod((NEW_LL_MATRIX[2][1] / C_0), 90.0));//(NEW_LL_MATRIX[2][1] / C_0);
            ROLL = atan(fmod((NEW_LL_MATRIX[2][0] / NEW_LL_MATRIX[2][2]), 180.0));//(NEW_LL_MATRIX[2][0] / NEW_LL_MATRIX[2][2]);
            YAW = atan(fmod((NEW_LL_MATRIX[0][1] / NEW_LL_MATRIX[1][1]), 360.0));//(NEW_LL_MATRIX[0][1] / NEW_LL_MATRIX[1][1]);

        }
        
        else{
            PITCH_0 = PITCH;
            ROLL_0 = ROLL;
            YAW_0 = YAW;
            //cout << wE << "    " << wN << "    " << wUp << endl;
            cout << PITCH_0 << "    " << ROLL_0 << "    " << YAW_0 << "    " << X << "    " << Y << endl;

            bodyToLocal(Acc_matrix_BL[0][0], Acc_matrix_BL[1][0], Acc_matrix_BL[2][0]);
            

            X = getX(Acc_matrix_ENUp[0][0]);
            Y = getY(Acc_matrix_ENUp[1][0]);
            //Z = getZ(Acc_matrix_ENUp[2][0]);

            //calculate angular velocity
            wE = -velocityY / R_latitude;
            wN = velocityX / R_longitude;
            wUp = velocityX / R_longitude * tan(latitude_0);

            Poisson();

            C_0 = sqrt(((NEW_LL_MATRIX[2][0] * NEW_LL_MATRIX[2][0]) + (NEW_LL_MATRIX[2][2] * NEW_LL_MATRIX[2][2])));

            PITCH = atan(fmod((NEW_LL_MATRIX[2][1] / C_0), 90.0));//(NEW_LL_MATRIX[2][1] / C_0);
            ROLL = atan(fmod((NEW_LL_MATRIX[2][0] / NEW_LL_MATRIX[2][2]), 180.0));//(NEW_LL_MATRIX[2][0] / NEW_LL_MATRIX[2][2]);
            YAW = atan(fmod((NEW_LL_MATRIX[0][1] / NEW_LL_MATRIX[1][1]), 360.0));//(NEW_LL_MATRIX[0][1] / NEW_LL_MATRIX[1][1]);
        }
    }

    Fin.close();

    cout << "Общее количество строк: " << numLines << endl;

    return 0;
}