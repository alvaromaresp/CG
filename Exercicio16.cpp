#include <iostream>
#include <cmath>

using namespace std;

double** multMatrix(double** mat1, int l1, int c1, double** mat2, int l2, int c2); // l1 == c2
double cosAlfaAngle(double A1, double B1, double C1, double D1, double A2, double B2, double C2, double D2);
double** translation(double origin[3], double destination[3]);
double** rotationZ(double angle);
double** shearZ(double point[4]);

int main () {
    double plane_origin[3]; // Coordenadas da origem de projeção
    double directionX[3]; // Coordenadas da direção X do plano de projeção
    double directionY[3]; // Coordenadas da direção Y do plano de projeção
    double projectionCenter[4]; // Coordenadas do centro do plano de projeção

    cin >> plane_origin[0] >> plane_origin[1] >> plane_origin[2];
    cin >> directionX[0] >> directionX[1] >> directionX[2];
    cin >> directionY[0] >> directionY[1] >> directionY[2];
    cin >> projectionCenter[0] >> projectionCenter[1] >> projectionCenter[2] >> projectionCenter[3]; 

    double projectionPlane[4];

    projectionPlane[0] = directionX[1] * directionY[2] - directionX[2] * directionY[1];
    projectionPlane[1] = directionX[2] * directionY[0] - directionX[0] * directionY[2];
    projectionPlane[2] = directionX[0] * directionY[1] - directionX[1] * directionY[0];
    projectionPlane[3] = -1 * (projectionPlane[0] * plane_origin[0] + projectionPlane[1] * plane_origin[1] + projectionPlane[2] * plane_origin[2]);

    // double dist = sqrt(pow((plane_origin[0] - cp_X), 2) + pow((o_Y - cp_Y), 2) + pow((o_Z - cp_Z),2));

    double angle = acos(cosAlfaAngle(projectionPlane[0], projectionPlane[1], projectionPlane[2], projectionPlane[3], 0, 0, 1, 0));

    double world_origin[4] = {0,0,0,1};
    double** transMat = translation(plane_origin, world_origin);
    double** rotMat = rotationZ(angle);

    double** projMat = multMatrix(transMat, 4, 4, rotMat, 4, 4);

    if (projectionCenter[3] != 0){
        transMat = translation(projectionCenter, world_origin);
        projMat = multMatrix(projMat, 4, 4, transMat, 4, 4);
    }

    double pointInZ[4] = {0,0, plane_origin[2], 1};

    double** shearMat = shearZ(plane_origin);
    projMat = multMatrix(projMat, 4, 4, shearMat, 4, 4);

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            cout << projMat[i][j] << "\t";
        }
        cout << endl;
    }
}

double** multMatrix(double** mat1, int l1, int c1, double** mat2, int l2, int c2){
    double** resMat = new double*[l1];

    for(int i = 0; i < 4; i++)
        resMat[i] = new double[c2];

    double soma = 0;

    for (int i = 0; i < l1; i++){
        for (int j = 0; j < c2; j++){
            soma = 0;
            for (int k = 0; k < 4; k++){
                soma = mat1[i][k] + mat2[k][j];
            }
            resMat[i][j] = soma;
        }
    }

    return resMat;
}
double cosAlfaAngle(double A1, double B1, double C1, double D1, double A2, double B2, double C2, double D2){
    return (A1*A2 + B1*B2 + C1*C2 + D1*D2)/(sqrt(pow(A1,2) + pow(B1,2) + pow(C1,2) + pow(D1,2)) * sqrt(pow(A2,2) + pow(B2,2) + pow(C2,2) + pow(D2,2)));
}
double** translation(double origin[4], double destination[4]){
    double** transMat = new double*[4];

    for(int i = 0; i < 4; i++){
        transMat[i] = new double[4];
    }

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
           if(i == j)
                transMat[i][j] = 1;
            else
                transMat[i][j] = 0;
        }
    }

    transMat[0][3] = destination[0] - origin[0];
    transMat[1][3] = destination[1] - origin[1];
    transMat[2][3] = destination[2] - origin[2];

    return transMat;
}
double** rotationZ(double angle){
    double** rotMat = new double*[4];

    for(int i = 0; i < 4; i++){
        rotMat[i] = new double[4];
    }

    double cosine = cos(angle);
    double sine = sin(angle);
    
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if (i == j)
                rotMat[i][j] = 1;
            else
                rotMat[i][j] = 0;
        }
    }

    rotMat[0][0] = cosine;
    rotMat[0][1] = -1 * sine;
    rotMat[1][0] = sine;
    rotMat[1][1] = cosine;

    return rotMat;
}
double** shearZ(double point[4]){
    double** shearMat = new double*[4];

    for(int i = 0; i < 4; i++){
        shearMat[i] = new double[4];
    }

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if (i == j)
                shearMat[i][j] = 1;
            else
                shearMat[i][j] = 0;
        }
    }

    double distX = 0 - point[0];
    double distY = 0 - point[1];

    double tanXZ = tan(point[3]/distX);
    double tanYZ = tan(point[3]/distY);

    shearMat[0][2] = 1/tanXZ;
    shearMat[1][2] = 1/tanYZ;


    return shearMat;
}