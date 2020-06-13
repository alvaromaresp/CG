#include <iostream>
#include <cmath>
#include <iomanip>      // std::setprecision

using namespace std;

double** multMatrix(double** mat1, int l1, int c1, double** mat2, int l2, int c2); // l1 == c2
double cosAlfaAngle(double A1, double B1, double C1, double A2, double B2, double C2);
double** translation(double origin[4], double destination[4]);
double** rotationY(double* xAxis);
double** rotationZ(double* xAxis);
double** rotationX(double* yAxis);
double** shear(double* point);
double scalarProduct(double* vector1, double* vector2);
void uniToBi(double** matrix, double* vector);

int main () {
    double* plane_origin = new double[4]; // Coordenadas da origem de projeção
    double* directionX = new double[4]; // Coordenadas da direção X do plano de projeção
    double* directionY = new double[4]; // Coordenadas da direção Y do plano de projeção
    double* projectionCenter = new double[4]; // Coordenadas do centro do plano de projeção

    // Input de dados
    cin >> plane_origin[0] >> plane_origin[1] >> plane_origin[2];
    cin >> directionX[0] >> directionX[1] >> directionX[2];
    cin >> directionY[0] >> directionY[1] >> directionY[2];
    cin >> projectionCenter[0] >> projectionCenter[1] >> projectionCenter[2] >> projectionCenter[3]; 

    // Coordenadas homogêneas
    plane_origin[3] = 1;
    directionX[3] = 0;
    directionY[3] = 0;

    // Vetor unidimensional para bidimensional
    double** dirXMat = new double*[4];
    for (int i = 0; i < 4; i++){
        dirXMat[i] = new double[1];
    }

    double** dirYMat = new double*[4];
    for (int i = 0; i < 4; i++){
        dirYMat[i] = new double[1];
    }

    double** projCenterMat = new double*[4];
    for (int i = 0; i < 4; i++){
        projCenterMat[i] = new double[1];
    }

    dirXMat[0][0] = directionX[0];
    dirXMat[1][0] = directionX[1];
    dirXMat[2][0] = directionX[2];
    dirXMat[3][0] = directionX[3];

    dirYMat[0][0] = directionY[0];
    dirYMat[1][0] = directionY[1];
    dirYMat[2][0] = directionY[2];
    dirYMat[3][0] = directionY[3];

    projCenterMat[0][0] = projectionCenter[0];
    projCenterMat[1][0] = projectionCenter[1];
    projCenterMat[2][0] = projectionCenter[2];
    projCenterMat[3][0] = projectionCenter[3];

    // Origem das coordenadas do mundo e translação para a mesma
    double world_origin[4] = {0,0,0,1};
    double** transMat = translation(plane_origin, world_origin);
    
    // Primeira rotação em Y levando o eixo X do plano pelo ângulo entre a projeção do mesmo em XZ e a direção X do mundo
    double** rotYMat = rotationY(directionX);

    dirXMat = multMatrix(rotYMat, 4, 4, dirXMat, 4, 1);
    dirYMat = multMatrix(rotYMat, 4, 4, dirYMat, 4, 1);
    projCenterMat = multMatrix(rotYMat, 4, 4, projCenterMat, 4, 1);

    cout << setprecision(5);
    cout << fixed;

    double** projMat = multMatrix(transMat, 4, 4, rotYMat, 4, 4);

    uniToBi(dirXMat, directionX);
    uniToBi(dirYMat, directionY);
    uniToBi(projCenterMat, projectionCenter);

    // Segunda rotação em Z levando a direção de X pelo ângulo entre a direção de X do plano e a direção X do mundo
    double** rotZMat = rotationZ(directionX);
    
    dirXMat = multMatrix(rotZMat, 4, 4, dirXMat, 4, 1);
    dirYMat = multMatrix(rotZMat, 4, 4, dirYMat, 4, 1);
    projCenterMat = multMatrix(rotZMat, 4, 4, projCenterMat, 4, 1);

    projMat = multMatrix(projMat, 4, 4, rotZMat, 4, 4);

    uniToBi(dirXMat, directionX);
    uniToBi(dirYMat, directionY);
    uniToBi(projCenterMat, projectionCenter);
    // Terceira rotação em Z levando a direção de X pelo ângulo entre a direção de X do plano e a direção X do mundo

    double** rotXMat = rotationX(directionY);

    dirXMat = multMatrix(rotXMat, 4, 4, dirXMat, 4, 1);
    dirYMat = multMatrix(rotXMat, 4, 4, dirYMat, 4, 1);
    projCenterMat = multMatrix(rotXMat, 4, 4, projCenterMat, 4, 1);

    projMat = multMatrix(projMat, 4, 4, rotXMat, 4, 4);
 
    double** shearMat;
    // Caso for projeção em perspectiva, transladar o CP para a origem do mundo e inclinar pelo centro
    if (projectionCenter[3] != 0){
        shearMat = shear(plane_origin);
        transMat = translation(projectionCenter, world_origin);
        projMat = multMatrix(projMat, 4, 4, transMat, 4, 4);
    } else {
        shearMat = shear(projectionCenter);
    }
    
    projMat = multMatrix(projMat, 4, 4, shearMat, 4, 4);
    cout << setw(8);
    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            cout << projMat[i][j] << " ";
        }
        cout << endl;
    }



}   

double** multMatrix(double** mat1, int r1, int c1, double** mat2, int r2, int c2){
    double** resMat = new double*[r1];

    for(int i = 0; i < r1; i++)
        resMat[i] = new double[c2];

    for(int i = 0; i < r1; i++)
        for (int j = 0; j < c2; j++)
            resMat[i][j] = 0.0;


    for (int i = 0; i < r1; i++){
        for (int j = 0; j < c2; j++){
            for (int k = 0; k < c1; k++){
                resMat[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return resMat;
}
double cosAlfaAngle(double A1, double B1, double C1, double A2, double B2, double C2){
    return ((A1*A2 + B1*B2 + C1*C2)/(sqrt(pow(A1,2) + pow(B1,2) + pow(C1,2)) * sqrt(pow(A2,2) + pow(B2,2) + pow(C2,2))));
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

    transMat[0][4] = destination[0] - origin[0];
    transMat[1][4] = destination[1] - origin[1];
    transMat[2][4] = destination[2] - origin[2];

    return transMat;
}
double** rotationY(double* xAxis){
    double** rotMat = new double*[4];

    for(int i = 0; i < 4; i++){
        rotMat[i] = new double[4];
    }

    double* xzNormal = new double[4];

    xzNormal[0] = 0;
    xzNormal[1] = 1;
    xzNormal[2] = 0;
    xzNormal[3] = 0;

    if (xzNormal[0] != 0)
        xzNormal[0] *= scalarProduct(xAxis, xzNormal) / ((pow(xzNormal[0], 2)) + (pow(xzNormal[1], 2)) + (pow(xzNormal[2], 2))); 
    if (xzNormal[1] != 0)
        xzNormal[1] *= scalarProduct(xAxis, xzNormal) / ((pow(xzNormal[0], 2)) + (pow(xzNormal[1], 2)) + (pow(xzNormal[2], 2)));  
    if (xzNormal[2] != 0)
        xzNormal[2] *= scalarProduct(xAxis, xzNormal) / ((pow(xzNormal[0], 2)) + (pow(xzNormal[1], 2)) + (pow(xzNormal[2], 2)));  
    
    double proj[3] = {0, 0, 0};

    proj[0] = xAxis[0] - xzNormal[0];
    proj[1] = xAxis[1] - xzNormal[1];
    proj[2] = xAxis[2] - xzNormal[2];


    double angle = acos(cosAlfaAngle(proj[0], proj[1], proj[2], 1, 0, 0));

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
    rotMat[0][2] = sine;
    rotMat[2][0] = -1  * sine;
    rotMat[2][2] = cosine;

    return rotMat;
}
double** rotationZ(double* xAxis){
    double** rotMat = new double*[4];

    for(int i = 0; i < 4; i++){
        rotMat[i] = new double[4];
    }

    double angle = acos(cosAlfaAngle(xAxis[0], xAxis[1], xAxis[2], 1, 0, 0));
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
    rotMat[0][1] = -1  * sine;
    rotMat[1][0] = sine;
    rotMat[1][1] = cosine;

    return rotMat;
}
double** rotationX(double* yAxis){
    double** rotMat = new double*[4];

    for(int i = 0; i < 4; i++){
        rotMat[i] = new double[4];
    }

    double angle = acos(cosAlfaAngle(yAxis[0], yAxis[1], yAxis[2], 0, 1, 0));
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
    rotMat[1][1] = cosine;
    rotMat[1][2] = -1  * sine;
    rotMat[2][1] = sine;
    rotMat[2][2] = cosine;

    return rotMat;
}
double** shear(double* point){
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

    if (point[2] != 0){
        shearMat[0][2] = (-1 * point[0]) / point[2];
        shearMat[1][2] = (-1 * point[1]) / point[2];
    }

    return shearMat;
}
double scalarProduct(double* vector1, double* vector2){
    return (vector1[0] * vector2[0]) + (vector1[1] * vector2[1]) + (vector1[2] * vector2[2]);
}
void uniToBi(double** matrix, double* vector){
    vector[0] = matrix [0][0];
    vector[1] = matrix [1][0];
    vector[2] = matrix [2][0];
    vector[3] = matrix [3][0];
}