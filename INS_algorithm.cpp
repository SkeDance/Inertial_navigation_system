#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

#define PI 3.14159265
#define ECC 0.081819190842622

using namespace std;

ofstream Fout("errors.txt");

bool alignment_flag = 0;
bool flag = 0;
bool flagX = 0;
bool flagY = 0;

const float g = 9.81523;
const double U = 15 * (PI / 180) / 3600; // Скорость вращения Земли, рад/сек
double dt = 0.005;
int t = 1;

// errors
double d_a = 1 * g * pow(10, -3);       // дрейф акселерометров
double w_dr = 1 * 0.1 / 3600;           // дрейф ДУСов

double w_M = 0.0002;                         // ошибка масштабных коэффициентов ДУСов
double a_M = 0.0005;                         // ошибка масштабных коэффициентов акселерометров

double err_Neort = 0.0003;//1.5 * pow(10, -3);   // ошибка неортогональности осей 
// 0,15мрад - test 
// 0,3мрад - calibrated
// 1,5мрад - uncalibrated

double FE;
double FN;
double FUp;
////

double CB_PL[3][3]; // матрица ошибок измерений
double CB_LL[3][3] = {0};

// shirota & dolgota
double fi_0 = 0;
double lambda_0 = 0;

double latitude = 0;
double longitude = 0;

double R_Earth = 6.371 * pow(10, 6);
double R_fi;
double R_lambda;

double C_0;

double wX, wY, wZ;
double aX, aY, aZ;

double PITCH_0 = 0;
double ROLL_0 = 0;
double YAW_0 = 0;

double PITCH = 0;
double ROLL = 0;
double YAW = 0;

// VN, North channel
double VN = 0;

double current_VN = 0;
double previous_VN = 0;
double integral_VN = 0;
double aNk = 0;

////

// VE, East channel
double VE = 0;

double current_VE = 0;
double previous_VE = 0;
double integral_VE = 0;
double aEk = 0;

////

// Vup, Up channel

double VUp = 0;
double aUPk = 0;

double wE = 0;
double wN = 0;
double wUp = 0;

double current_aX = 0;
double previous_aX = 0;

double current_aY = 0;
double previous_aY = 0;

double current_aZ = 0;
double previous_aZ = 0;

double velocityX = 0;
double prev_velocityX = 0;

double velocityY = 0;
double prev_velocityY = 0;

double velocityZ = 0;
double prev_velocityZ = 0;

double fi, lambda, H;

double matrix_LL[3][3]; // матрица для выставки на первом такте

double Acc_matrix_BL[3][1];  // матрица с показаниями акселерометров
double Gyro_matrix_BL[3][1]; // матрица с показаниями ДУСов

double Acc_matrix_ENUp[3][1]; // матрица ускорений в системе ENUp

double matrix_W_B[3][3];  // 3x3 матрица показаний ДУСов
double matrix_W_LL[3][3]; // матрица посчитанных угловых скоростей из показаний акслерометра

double matrix_first[3][3];  // первый элемент уравнения Пуассона
double matrix_second[3][3]; // второй элемент уравнения Пуассона

double NEW_LL_MATRIX[3][3]; // новая матрица, которая используется для перевода в систему ENUp на последующих тактах работы алгоритма

double Gyro_matrix_BL_E[3][1]; // матрица показаний ДУСов с учетом ошибок
double Acc_matrix_BL_E[3][1];  // матрица показаний акселерометров с учетом ошибок

double ErrorsMatrixDUS[3][3] = {
    {1.0 + w_M, err_Neort, err_Neort},
    {err_Neort, 1.0 + w_M, err_Neort},
    {err_Neort, err_Neort, 1.0 + w_M}
};

double ErrorsMatrixAcc[3][3] = {
    {1.0 + a_M, err_Neort, err_Neort},
    {err_Neort, 1.0 + a_M, err_Neort},
    {err_Neort, err_Neort, 1.0 + a_M}
};

float takt = 0;

void alignment()
{
    double g_0 = sqrt(pow(Acc_matrix_BL[0][0], 2) + pow(Acc_matrix_BL[1][0], 2) + pow(Acc_matrix_BL[2][0], 2));
    double w_dzetta = U * sin(fi_0);
    double w_nu = U * cos(fi_0);
    matrix_LL[2][0] = -Acc_matrix_BL[0][0] / g_0;
    matrix_LL[2][1] = -Acc_matrix_BL[1][0] / g_0;
    matrix_LL[2][2] = -Acc_matrix_BL[2][0] / g_0;
    matrix_LL[1][0] = (Gyro_matrix_BL[0][0] - w_dzetta * matrix_LL[2][0]);
    matrix_LL[1][1] = (Gyro_matrix_BL[1][0] - w_dzetta * matrix_LL[2][1]);
    matrix_LL[1][2] = (Gyro_matrix_BL[2][0] - w_dzetta * matrix_LL[2][2]);
    matrix_LL[0][0] = (matrix_LL[1][1] * matrix_LL[2][2] - matrix_LL[1][2] * matrix_LL[2][1]);
    matrix_LL[0][1] = (matrix_LL[1][2] * matrix_LL[2][0] - matrix_LL[1][0] * matrix_LL[2][2]);
    matrix_LL[0][2] = (matrix_LL[1][0] * matrix_LL[2][1] - matrix_LL[1][1] * matrix_LL[2][0]);
}

double normalizeAngle(double degrees)
{
    const double rounding_threshold = 1e-2; // Порог 0.01°

    while (degrees < 0.0)
    {
        degrees += 360.0;
    }

    while (degrees >= 360.0)
    {
        degrees -= 360.0;
    }

    if (360.0 - degrees < rounding_threshold)
    {
        return 0.0;
    }

    return degrees;
}

double DegreesToRads(double degree)
{
    return degree * PI / 180;
}

double RadsToDegrees(double rads)
{
    return rads * 180 / PI;
}

double getSpeedVE(double acceleration)
{
    previous_aX = current_aX;
    current_aX = acceleration;
    VE += ((((current_aX + previous_aX) / 2.0) + aEk) * dt);
    return VE;
}

double getSpeedVN(double acceleration)
{
    previous_aY = current_aY;
    current_aY = acceleration;
    VN += ((((current_aY + previous_aY) / 2.0) + aNk) * dt);
    return VN;
}

double getSpeedVUp(double acceleration)
{
    previous_aZ = current_aZ;
    current_aZ = acceleration;
    VUp += ((((current_aZ + previous_aZ) / 2.0) + aUPk) * dt); //- wE * VN + wN * VE + U * cos(fi_0) * VE - g) * dt;
    return VUp;
}

double getFi(double speedVN)
{
    previous_VN = current_VN;
    current_VN = speedVN;
    double avg_speed = ((current_VN + previous_VN) / 2.0);
    integral_VN += ((avg_speed / (R_fi + H)) * dt);
    return fi_0 + integral_VN;
}

double getLambda(double speedVE)
{
    previous_VE = current_VE;
    current_VE = speedVE;
    double avg_speed = ((current_VE + previous_VE) / 2.0);
    integral_VE += ((avg_speed / ((R_lambda + H) * cos(fi))) * dt);
    return lambda_0 + integral_VE;
}

void UpdateLLMatrix(double LL[3][3], double LL_dt[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            LL[i][j] += LL_dt[i][j] * dt;
        }
    }
}

void normalize(double C_B_LL[3][3], int t)
{
    for (int i = 0; i < 3; i++)
    {
        double s;

        if (t % 2 == 0)
        { // Если такт четный, нормализация по строкам
            s = sqrt(C_B_LL[i][0] * C_B_LL[i][0] +
                     C_B_LL[i][1] * C_B_LL[i][1] +
                     C_B_LL[i][2] * C_B_LL[i][2]);

            C_B_LL[i][0] /= s;
            C_B_LL[i][1] /= s;
            C_B_LL[i][2] /= s;
        }
        else
        { // Если такт нечетный, нормализация по столбцам
            s = sqrt(C_B_LL[0][i] * C_B_LL[0][i] +
                     C_B_LL[1][i] * C_B_LL[1][i] +
                     C_B_LL[2][i] * C_B_LL[2][i]);

            C_B_LL[0][i] /= s;
            C_B_LL[1][i] /= s;
            C_B_LL[2][i] /= s;
        }
    }
}

void Poisson(const double a[3][3], const double b[3][3], double result[3][3])
{
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            result[i][j] = 0; // Инициализация элемента результата
            for (int k = 0; k < 3; ++k)
            {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void addErrorsDUS(const double A[3][3], const double B[3][1], double C[3][1]) {
    for (int i = 0; i < 3; i++) {
        C[i][0] = 0;
        for (int k = 0; k < 3; k++) {
            C[i][0] += B[k][0] * A[i][k]; // Суммируем все компоненты
        }
        C[i][0] += DegreesToRads(w_dr); // Добавляем дрейф после суммирования
    }
}

void addErrorsAcc(const double A[3][3], const double B[3][1], double C[3][1]) {
    for (int i = 0; i < 3; i++) {
        C[i][0] = 0;
        for (int k = 0; k < 3; k++) {
            C[i][0] += B[k][0] * A[i][k]; // Суммируем все компоненты
        }
        C[i][0] += d_a; // Добавляем дрейф
    }
}

void bodyToLocal(const double A[3][3], const double B[3][1], double C[3][1])
{
    for (int i = 0; i < 3; i++)
    {
        C[i][0] = 0; // Инициализация элемента результата
        for (int k = 0; k < 3; k++)
        {
            C[i][0] += A[i][k] * B[k][0];
        }
    }
}

void subtractMatrices(const double matrix1[3][3], const double matrix2[3][3], double result[3][3])
{
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            result[i][j] = matrix1[i][j] - matrix2[i][j];
        }
    }
}

void matrix()
{
    matrix_LL[0][0] = cos(DegreesToRads(ROLL_0)) * cos(DegreesToRads(YAW_0)) + sin(DegreesToRads(PITCH_0)) * sin(DegreesToRads(ROLL_0)) * sin(DegreesToRads(YAW_0));
    matrix_LL[0][1] = -cos(DegreesToRads(ROLL_0)) * sin(DegreesToRads(YAW_0)) + sin(DegreesToRads(PITCH_0)) * sin(DegreesToRads(ROLL_0)) * cos(DegreesToRads(YAW_0));
    matrix_LL[0][2] = cos(DegreesToRads(PITCH_0)) * sin(DegreesToRads(ROLL_0));
    matrix_LL[1][0] = cos(DegreesToRads(PITCH_0)) * sin(DegreesToRads(YAW_0));
    matrix_LL[1][1] = cos(DegreesToRads(PITCH_0)) * cos(DegreesToRads(YAW_0));
    matrix_LL[1][2] = sin(DegreesToRads(PITCH_0));
    matrix_LL[2][0] = sin(DegreesToRads(ROLL_0)) * cos(DegreesToRads(YAW_0)) - sin(DegreesToRads(PITCH_0)) * cos(DegreesToRads(ROLL_0)) * sin(DegreesToRads(YAW_0));
    matrix_LL[2][1] = -sin(DegreesToRads(ROLL_0)) * sin(DegreesToRads(YAW_0)) - sin(DegreesToRads(PITCH_0)) * cos(DegreesToRads(ROLL_0)) * sin(DegreesToRads(YAW_0));
    matrix_LL[2][2] = cos(DegreesToRads(PITCH_0)) * cos(DegreesToRads(ROLL_0));
}

void matrix_W_DUS()
{
    matrix_W_B[0][0] = 0;
    matrix_W_B[0][1] = -Gyro_matrix_BL_E[2][0]; // rad/sec
    matrix_W_B[0][2] = Gyro_matrix_BL_E[1][0];  // rad/sec
    matrix_W_B[1][0] = Gyro_matrix_BL_E[2][0];  // rad/sec
    matrix_W_B[1][1] = 0;
    matrix_W_B[1][2] = -Gyro_matrix_BL_E[0][0]; // rad/sec
    matrix_W_B[2][0] = -Gyro_matrix_BL_E[1][0]; // rad/sec
    matrix_W_B[2][1] = Gyro_matrix_BL_E[0][0];  // rad/sec
    matrix_W_B[2][2] = 0;
}

void matrix_W_ENUp()
{
    matrix_W_LL[0][0] = 0;
    matrix_W_LL[0][1] = -wUp;
    matrix_W_LL[0][2] = wN + (0.0024 / 360 * PI / 180);
    matrix_W_LL[1][0] = wUp;
    matrix_W_LL[1][1] = 0;
    matrix_W_LL[1][2] = -(wE + (-0.0001 / 360 * PI / 180));
    matrix_W_LL[2][0] = -(wN + (0.0024 / 360 * PI / 180));
    matrix_W_LL[2][1] = (wE + (-0.0001 / 360 * PI / 180));
    matrix_W_LL[2][2] = 0;
}

void set()
{
    matrix();
    FE = d_a / g;
    FN = d_a / g;
    FUp = w_dr / (U * cos(fi_0));

    CB_PL[0][0] = 1;
    CB_PL[0][1] = -FUp;
    CB_PL[0][2] = FN;
    CB_PL[1][0] = FUp;
    CB_PL[1][1] = 1;
    CB_PL[1][2] = -FE;
    CB_PL[2][0] = -FN;
    CB_PL[2][1] = FE;
    CB_PL[2][2] = 1;

    Poisson(matrix_LL, CB_PL, CB_LL);
}

int main()
{
    string path = "20170922_Cessna172_200Hz_Ref.txt";
    ifstream Fin;
    Fin.open(path);

    if (!Fout.is_open())
    {
        cerr << "Ошибка создания файла для записи!" << endl;
        return 1;
    }

    // // Формат данных в файле записи
    Fout << "Time\tPITCH\tROLL\tYAW\tLatitude\tLongitude\n";

    if (!Fin.is_open())
    {
        cerr << "Ошибка открытия файла!" << endl;
        return 1;
    }

    string str; // Переменная для записи строк
    int numLines = 0;

    while (getline(Fin, str))
    {
        ++numLines;

        // Разбиение на строки при помощи std::istringstream
        istringstream iss(str);
        string token;
        int tokenCount = 0;

        while (iss >> token)
        {
            ++tokenCount;
            if (tokenCount == 1)
            {
                // time
                takt = stof(token);
            }
            if (tokenCount == 2)
            {
                // wX
                Gyro_matrix_BL[0][0] = DegreesToRads(stof(token));// + DegreesToRads(w_dr); // rad/sec
            }

            if (tokenCount == 3)
            {
                // wY
                Gyro_matrix_BL[1][0] = DegreesToRads(stof(token));// + DegreesToRads(w_dr); // rad/sec
            }

            if (tokenCount == 4)
            {
                // wZ
                Gyro_matrix_BL[2][0] = DegreesToRads(stof(token));// + DegreesToRads(w_dr); // rad/sec
            }

            if (tokenCount == 5)
            {
                // aX
                Acc_matrix_BL[0][0] = stof(token);// + d_a;
            }

            if (tokenCount == 6)
            {
                // aY
                Acc_matrix_BL[1][0] = stof(token);// + d_a;
            }

            if (tokenCount == 7)
            {
                // aZ
                Acc_matrix_BL[2][0] = stof(token);// + d_a;
            }

            if (tokenCount == 15)
            {
                // Height
                H = stof(token); // metres
            }

            if (alignment_flag == 0)
            {
                if (tokenCount == 8)
                {
                    // крен
                    ROLL_0 = stof(token); // degrees
                }

                if (tokenCount == 9)
                {
                    // тангаж
                    PITCH_0 = stof(token); // degrees
                }

                if (tokenCount == 10)
                {
                    // курс
                    YAW_0 = stof(token); // degrees
                }

                if (tokenCount == 13)
                {
                    // широта
                    fi_0 = DegreesToRads(stof(token)); // rads
                }

                if (tokenCount == 14)
                {
                    // долгота
                    lambda_0 = DegreesToRads(stof(token)); // rads
                }
            }
        }

        // alignment
        if (t < 40296) // t < 48060 |||| (takt <  201.480000) == (t < 40296)
        {
            if (alignment_flag == 0)
            {
                // Добавление ошибок 
                addErrorsDUS(ErrorsMatrixDUS, Gyro_matrix_BL, Gyro_matrix_BL_E);
                addErrorsAcc(ErrorsMatrixAcc, Acc_matrix_BL, Acc_matrix_BL_E);

                set();
                // matrix();
                bodyToLocal(matrix_LL, Acc_matrix_BL_E, Acc_matrix_ENUp);

                VE = getSpeedVE(Acc_matrix_ENUp[0][0]);
                VN = getSpeedVN(Acc_matrix_ENUp[1][0]);
                VUp = getSpeedVUp(Acc_matrix_ENUp[2][0]);

                aEk = ((-wN * VUp) + (wUp * VN) +(-U * cos(fi_0) * VUp) + (U * sin(fi_0) * VN));
                aNk = ((-wE * VUp) + (-wUp * VE) + (-U * sin(fi_0) * VE));
                aUPk = ((-wE * VN) + (wN * VE) + (U * cos(fi_0) * VE - g));

                R_fi = (R_Earth * (1 - pow(ECC, 2))) / pow(1 - pow(ECC, 2) * pow(sin(fi_0), 2), 3.0 / 2.0);
                R_lambda = (R_Earth) / pow(1 - pow(ECC, 2) * pow(sin(fi_0), 2), 1.0 / 2.0);

                // calculate angular velocity
                wE = -VN / (R_fi + H);                                 // rad/s
                wN = VE / (R_lambda + H) + U * cos(fi_0);              // rad/s
                wUp = VE / (R_lambda + H) * tan(fi_0) + U * sin(fi_0); // rad/s

                fi = getFi(VN);         // rads
                lambda = getLambda(VE); // rads

                matrix_W_ENUp(); // составляем матрицу рассчитанных угловых скоростей
                matrix_W_DUS();  // составляем матрицу угловых скоростей из показаний ДУСов

                Poisson(matrix_LL, matrix_W_B, matrix_first);                 // расчет - первый элемент уравнения Пуассона
                Poisson(matrix_W_LL, matrix_LL, matrix_second);               // расчет - второй элемент уравнеия Пуассона
                subtractMatrices(matrix_first, matrix_second, NEW_LL_MATRIX); // решение уравнения Пуассона

                UpdateLLMatrix(matrix_LL, NEW_LL_MATRIX);

                normalize(matrix_LL, t);
                t++; // такт работы

                // Ориентация
                C_0 = sqrt(matrix_LL[2][0] * matrix_LL[2][0] + matrix_LL[2][2] * matrix_LL[2][2]);
                PITCH = RadsToDegrees(atan(matrix_LL[2][1] / C_0));
                ROLL = RadsToDegrees(-atan2(matrix_LL[2][0], matrix_LL[2][2]));
                YAW = normalizeAngle(RadsToDegrees(atan2(matrix_LL[0][1], matrix_LL[1][1])));

                std::cout << "широта   " << RadsToDegrees(fi) << "  долгота    " << RadsToDegrees(lambda) << "  крен  " << PITCH << "  тангаж   " << ROLL << "  курс  " << YAW << endl;

                Fout << takt << "\t"
                    << PITCH << "\t"
                    << ROLL << "\t"
                    << YAW << "\t"
                    << RadsToDegrees(fi) << "\t"
                    << RadsToDegrees(lambda) << "\n";

                alignment_flag = 1;
                
            }

            else
            {   
                // Добавление ошибок 
                addErrorsDUS(ErrorsMatrixDUS, Gyro_matrix_BL, Gyro_matrix_BL_E);
                addErrorsAcc(ErrorsMatrixAcc, Acc_matrix_BL, Acc_matrix_BL_E);

                bodyToLocal(matrix_LL, Acc_matrix_BL_E, Acc_matrix_ENUp);

                VE = getSpeedVE(Acc_matrix_ENUp[0][0]);
                VN = getSpeedVN(Acc_matrix_ENUp[1][0]);
                VUp = getSpeedVUp(Acc_matrix_ENUp[2][0]);

                aEk = ((-wN * VUp) + (wUp * VN) +(-U * cos(fi) * VUp) + (U * sin(fi) * VN));
                aNk = ((-wE * VUp) + (-wUp * VE) + (-U * sin(fi) * VE));
                aUPk = ((-wE * VN) + (wN * VE) + (U * cos(fi) * VE - g));

                R_fi = (R_Earth * (1 - pow(ECC, 2))) / pow(1 - pow(ECC, 2) * pow(sin(fi), 2), 3.0 / 2.0);
                R_lambda = (R_Earth) / pow(1 - pow(ECC, 2) * pow(sin(fi_0), 2), 1.0 / 2.0);

                // calculate angular velocity
                wE = -VN / (R_fi + H);                                 // rad/s
                wN = VE / (R_lambda + H) + U * cos(fi);              // rad/s
                wUp = VE / (R_lambda + H) * tan(fi) + U * sin(fi); // rad/s

                fi = getFi(VN);         // rads
                lambda = getLambda(VE); // rads

                matrix_W_ENUp(); // составляем матрицу рассчитанных угловых скоростей
                matrix_W_DUS();  // составляем матрицу угловых скоростей из показаний ДУСов

                Poisson(matrix_LL, matrix_W_B, matrix_first);                 // расчет - первый элемент уравнения Пуассона
                Poisson(matrix_W_LL, matrix_LL, matrix_second);               // расчет - второй элемент уравнеия Пуассона
                subtractMatrices(matrix_first, matrix_second, NEW_LL_MATRIX); // решение уравнения Пуассона

                UpdateLLMatrix(matrix_LL, NEW_LL_MATRIX);

                normalize(matrix_LL, t);
                t++; // такт работы

                // Ориентация
                C_0 = sqrt(matrix_LL[2][0] * matrix_LL[2][0] + matrix_LL[2][2] * matrix_LL[2][2]);
                PITCH = RadsToDegrees(atan(matrix_LL[2][1] / C_0));
                ROLL = RadsToDegrees(-atan2(matrix_LL[2][0], matrix_LL[2][2]));
                YAW = normalizeAngle(RadsToDegrees(atan2(matrix_LL[0][1], matrix_LL[1][1])));

                std::cout << "широта   " << RadsToDegrees(fi) << "  долгота    " << RadsToDegrees(lambda) << "  крен  " << PITCH << "  тангаж   " << ROLL << "  курс  " << YAW << endl;

                Fout << takt << "\t"
                    << PITCH << "\t"
                    << ROLL << "\t"
                    << YAW << "\t"
                    << RadsToDegrees(fi) << "\t"
                    << RadsToDegrees(lambda) << "\n";

            }
        }

        else
        {   
            // Добавление ошибок 
            addErrorsDUS(ErrorsMatrixDUS, Gyro_matrix_BL, Gyro_matrix_BL_E);
            addErrorsAcc(ErrorsMatrixAcc, Acc_matrix_BL, Acc_matrix_BL_E);

            bodyToLocal(matrix_LL, Acc_matrix_BL_E, Acc_matrix_ENUp);

            VE = getSpeedVE(Acc_matrix_ENUp[0][0]);
            VN = getSpeedVN(Acc_matrix_ENUp[1][0]);
            VUp = getSpeedVUp(Acc_matrix_ENUp[2][0]);

            aEk = ((-wN * VUp) + (wUp * VN) +(-U * cos(fi) * VUp) + (U * sin(fi) * VN));
            aNk = ((-wE * VUp) + (-wUp * VE) + (-U * sin(fi) * VE));
            aUPk = ((-wE * VN) + (wN * VE) + (U * cos(fi) * VE - g));

            // calculate angular velocity
            wE = -VN / (R_fi + H);                             // rad/s
            wN = VE / (R_lambda + H) + U * cos(fi);            // rad/s
            wUp = VE / (R_lambda + H) * tan(fi) + U * sin(fi); // rad/s

            // Расчет радиусов
            R_fi = (R_Earth * (1 - pow(ECC, 2))) / pow(1 - pow(ECC, 2) * pow(sin(fi), 2), 3.0 / 2.0);
            R_lambda = (R_Earth) / pow(1 - pow(ECC, 2) * pow(sin(fi), 2), 1.0 / 2.0);

            // Навигация
            fi = getFi(VN);         // rads
            lambda = getLambda(VE); // rads

            matrix_W_ENUp(); // составляем матрицу рассчитанных угловых скоростей
            matrix_W_DUS();  // составляем матрицу угловых скоростей из показаний ДУСов

            Poisson(matrix_LL, matrix_W_B, matrix_first);                 // расчет - первый элемент уравнения Пуассона
            Poisson(matrix_W_LL, matrix_LL, matrix_second);               // расчет - второй элемент уравнеия Пуассона
            subtractMatrices(matrix_first, matrix_second, NEW_LL_MATRIX); // решение уравнения Пуассона

            UpdateLLMatrix(matrix_LL, NEW_LL_MATRIX);

            normalize(matrix_LL, t);
            t++; // такт работы

            // Ориентация
            C_0 = sqrt(matrix_LL[2][0] * matrix_LL[2][0] + matrix_LL[2][2] * matrix_LL[2][2]);
            PITCH = RadsToDegrees(atan(matrix_LL[2][1] / C_0));
            ROLL = RadsToDegrees(-atan2(matrix_LL[2][0], matrix_LL[2][2]));
            YAW = normalizeAngle(RadsToDegrees(atan2(matrix_LL[0][1], matrix_LL[1][1])));

            std::cout << "широта   " << RadsToDegrees(fi) << "  долгота    " << RadsToDegrees(lambda) << "  крен  " << PITCH << "  тангаж   " << ROLL << "  курс  " << YAW << endl;

            Fout << takt << "\t"
                 << PITCH << "\t"
                 << ROLL << "\t"
                 << YAW << "\t"
                 << RadsToDegrees(fi) << "\t"
                 << RadsToDegrees(lambda) << "\n";
        }
    }

    Fin.close();  // Закрываем файл после чтения
    Fout.close(); // Закрываем файл после записи

    std::cout << "Общее количество строк: " << numLines << endl;

    return 0;
}