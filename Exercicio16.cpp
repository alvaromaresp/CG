#include <iostream>
#include <cmath>
#include <iomanip>      // std::setprecision

using namespace std;


double** multMatrix(double** mat1, int l1, int c1, double** mat2, int l2, int c2);
double cosAlfaAngle(double A1, double B1, double C1, double A2, double B2, double C2);
double** translation(double origin[4], double destination[4]);
double** rotationY(double** xAxis);
double** rotationZ(double** xAxis);
double** rotationX(double** yAxis);
double** shearToZ(double** directionX, double** directionY);
double scalarProduct(double** vector1, double** vector2);
void uniToBi(double** matrix, double* vector);
void initializeMatrix(double** matrix, int n, int m);
void perspectiveProjection(double** directionX, double** directionY, double** projectionCenter, double** planeOriginMat, double** world_origin);
void paralalelProjection(double** directionX, double** directionY, double** projectionVector, double** planeOriginMat, double** world_origin);
void translateValues(double** transMat, double** directionX, double** directionY, double** projectionCenter, double** planeOriginMat);
void updateValues(double** transMat, double** directionX, double** directionY, double** projectionCenter);
void printMatrix(double** matrix, int n, int m);

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

    double** planeOriginMat = new double*[4];
    for (int i = 0; i < 4; i++){
        planeOriginMat[i] = new double[1];
    }


    double** world_origin = new double*[4];
    for (int i = 0; i < 4; i++){
        world_origin[i] = new double[1];
    }


    initializeMatrix(world_origin, 4, 1);
    world_origin[3][0] = 1;


    initializeMatrix(dirXMat, 4, 1);
    initializeMatrix(dirYMat, 4, 1);
    initializeMatrix(projCenterMat, 4, 1);
    initializeMatrix(planeOriginMat, 4, 1);
    

    uniToBi(dirXMat, directionX);
    uniToBi(dirYMat, directionY);
    uniToBi(projCenterMat, projectionCenter);
    uniToBi(planeOriginMat, plane_origin);

    delete[] directionX;
    delete[] directionY;
    delete[] projectionCenter;
    delete[] plane_origin;


    if (projCenterMat[3][0] == 0)
        paralalelProjection(dirXMat, dirYMat, projCenterMat, planeOriginMat, world_origin);
    else
        perspectiveProjection(dirXMat, dirYMat, projCenterMat, planeOriginMat, world_origin);

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
double** translation(double** origin, double** destination){

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

    transMat[0][3] = destination[0][0] - origin[0][0];
    transMat[1][3] = destination[1][0] - origin[1][0];
    transMat[2][3] = destination[2][0] - origin[2][0];

    return transMat;
}
double** rotationY(double** xAxis){
    double** rotMat = new double*[4];

    for(int i = 0; i < 4; i++){
        rotMat[i] = new double[4];
    }

    double** xzNormal = new double*[4];

    for(int i = 0; i < 4; i++){
        xzNormal[i] = new double[1];
    }

    xzNormal[0][0] = 0;
    xzNormal[1][0] = 1;
    xzNormal[2][0] = 0;
    xzNormal[3][0] = 0;

    if (xzNormal[0][0] != 0)
        xzNormal[0][0] *= scalarProduct(xAxis, xzNormal) / ((pow(xzNormal[0][0], 2)) + (pow(xzNormal[1][0], 2)) + (pow(xzNormal[2][0], 2))); 
    if (xzNormal[1][0] != 0)
        xzNormal[1][0] *= scalarProduct(xAxis, xzNormal) / ((pow(xzNormal[0][0], 2)) + (pow(xzNormal[1][0], 2)) + (pow(xzNormal[2][0], 2)));  
    if (xzNormal[2][0] != 0)
        xzNormal[2][0] *= scalarProduct(xAxis, xzNormal) / ((pow(xzNormal[0][0], 2)) + (pow(xzNormal[1][0], 2)) + (pow(xzNormal[2][0], 2)));  
    
    double proj[3] = {0,0,0};

    proj[0] = xAxis[0][0] - xzNormal[0][0];
    proj[1] = xAxis[1][0] - xzNormal[1][0];
    proj[2] = xAxis[2][0] - xzNormal[2][0];


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
double** rotationZ(double** xAxis){
    double** rotMat = new double*[4];

    for(int i = 0; i < 4; i++){
        rotMat[i] = new double[4];
    }

    double angle = acos(cosAlfaAngle(xAxis[0][0], xAxis[1][0], xAxis[2][0], 1, 0, 0));
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
double** rotationX(double** yAxis){
    double** rotMat = new double*[4];

    for(int i = 0; i < 4; i++){
        rotMat[i] = new double[4];
    }

    double angle = acos(cosAlfaAngle(yAxis[0][0], yAxis[1][0], yAxis[2][0], 0, 1, 0));
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
double** shearToZ(double** directionX, double** directionY){
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

    double x = (directionX[1][0] * directionY[2][0]) - (directionX[2][0] * directionY[1][0]);
    double y = (directionX[2][0] * directionY[0][0]) - (directionX[0][0] * directionY[2][0]);
    double z = (directionX[0][0] * directionY[1][0]) - (directionX[1][0] * directionY[0][0]);

    if (z != 0){
        shearMat[0][2] = -1 * x / z;
        shearMat[1][2] = -1 * y / z;
    }

    return shearMat;
}
double** shearParalelToZ(double** directionX, double** directionY){
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

    double x = (directionX[1][0] * directionY[2][0]) - (directionX[2][0] * directionY[1][0]);
    double y = (directionX[2][0] * directionY[0][0]) - (directionX[0][0] * directionY[2][0]);
    double z = (directionX[0][0] * directionY[1][0]) - (directionX[1][0] * directionY[0][0]);

    double angle = acos(cos(cosAlfaAngle(x, y, z, 0, 0, 1)));

    if (angle != 0){
        shearMat[0][2] = 1 / tan(angle);
        shearMat[1][2] = 1 / tan(angle);
    }

    return shearMat;
}
double scalarProduct(double** vector1, double** vector2){
    return (vector1[0][0] * vector2[0][0]) + (vector1[1][0] * vector2[1][0]) + (vector1[2][0] * vector2[2][0]);
}
void uniToBi(double** matrix, double* vector){
    matrix [0][0] = vector[0];
    matrix [1][0] = vector[1];
    matrix [2][0] = vector[2];
    matrix [3][0] = vector[3];
}
void initializeMatrix(double** matrix, int n, int m){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
                matrix[i][j] = 0;
        }
    }
}

void perspectiveProjection(double** directionX, double** directionY, double** projectionCenter, double** planeOriginMat, double** world_origin){
    double** projMat = NULL;
    double** rotYMat = NULL;
    double** rotZMat = NULL;
    double** rotXMat = NULL;
    double** transMat = NULL;
    double** shearMat = NULL;;

    // 1: Transladar a origem do PP para a origem do SCM

    projMat = translation(planeOriginMat, world_origin);
    translateValues(projMat, directionX, directionY, projectionCenter, planeOriginMat);

    // 2: Rotacionar o plano de projeção

    // 2.1: Rotação em Y
    rotYMat = rotationY(directionX);
    projMat = multMatrix(projMat, 4, 4, rotYMat, 4, 4);
    updateValues(rotYMat, directionX, directionY, projectionCenter);

    // 2.2: Rotação em Z
    rotZMat = rotationZ(directionX);
    projMat = multMatrix(projMat, 4, 4, rotZMat, 4, 4);
    updateValues(rotZMat, directionX, directionY, projectionCenter);

    // 2.3: Rotação em X
    rotXMat = rotationX(directionY);
    projMat = multMatrix(projMat, 4, 4, rotXMat, 4, 4);
    updateValues(rotXMat, directionX, directionY, projectionCenter);

    // 3: Transladar CP para origem
    transMat = translation(projectionCenter, world_origin);
    projMat = multMatrix(projMat, 4, 4, transMat, 4, 4);
    translateValues(transMat, directionX, directionY, projectionCenter, planeOriginMat);

    // 4: Shear até que a origem do plano de projeção esteja no eixo Z
    shearMat = shearToZ(directionX, directionY);
    projMat = multMatrix(projMat, 4, 4, shearMat, 4, 4);
    translateValues(shearMat, directionX, directionY, projectionCenter, planeOriginMat);

    printMatrix(projMat, 4, 4);

    delete[] projMat;
    delete[] rotYMat;
    delete[] rotXMat;
    delete[] rotZMat;
    delete[] shearMat;
    delete[] transMat;
}

void paralalelProjection(double** directionX, double** directionY, double** projectionVector, double** planeOriginMat, double** world_origin){
    double** projMat = NULL;
    double** rotYMat = NULL;
    double** rotZMat = NULL;
    double** rotXMat = NULL;
    double** shearMat = NULL;;
    // cout << "Projeção paralela" << endl;
    // 1: Transladar a OP para a origem do SCM
    
    projMat = translation(planeOriginMat, world_origin);
    translateValues(projMat, directionX, directionY, projectionVector, planeOriginMat);

    // cout << " Etapa 1 concluída" << endl;
    // 2: Rotacionar o PP

    // 2.1: Rotação em Y

    rotYMat = rotationY(directionX);
    projMat = multMatrix(projMat, 4, 4, rotYMat, 4, 4);
    updateValues(rotYMat, directionX, directionY, projectionVector);

    // printMatrix(projMat, 4, 4);

    // cout << " Etapa 2.1 concluída" << endl;
    // 2.2: Rotação em Z

    rotZMat = rotationZ(directionX);
    projMat = multMatrix(projMat, 4, 4, rotZMat, 4, 4);
    updateValues(rotZMat, directionX, directionY, projectionVector);

    // printMatrix(projMat, 4, 4);

    // cout << " Etapa 2.2 concluída" << endl;
    // 2.3: Rotação em X

    rotXMat = rotationX(directionY);
    projMat = multMatrix(projMat, 4, 4, rotXMat, 4, 4);
    updateValues(rotXMat, directionX, directionY, projectionVector);

    // printMatrix(projMat, 4, 4);

    // cout << " Etapa 2.3 concluída" << endl;

    // 3: Shear até que o vetor de projeção esteja paralelo ao eixo Z
    shearMat = shearParalelToZ(directionX, directionY);
    projMat = multMatrix(projMat, 4, 4, shearMat, 4, 4);
    translateValues(shearMat, directionX, directionY, projectionVector, planeOriginMat);

    // printMatrix(projMat, 4, 4);

    // cout << " Etapa 3 concluída" << endl;

    printMatrix(projMat, 4, 4);

    delete[] projMat;
    delete[] rotYMat;
    delete[] rotXMat;
    delete[] rotZMat;
    delete[] shearMat;
}
void translateValues(double** transMat, double** directionX, double** directionY, double** projectionCenter, double** planeOriginMat){
    projectionCenter = multMatrix(transMat, 4, 4, projectionCenter, 4, 1);
    planeOriginMat = multMatrix(transMat, 4, 4, planeOriginMat, 4, 1);
    directionX = multMatrix(transMat, 4, 4, directionX, 4, 1);
    directionY = multMatrix(transMat, 4, 4, directionY, 4, 1);
}
void updateValues(double** transMat, double** directionX, double** directionY, double** projectionCenter){
    projectionCenter = multMatrix(transMat, 4, 4, projectionCenter, 4, 1);
    directionX = multMatrix(transMat, 4, 4, directionX, 4, 1);
    directionY = multMatrix(transMat, 4, 4, directionY, 4, 1);
}

void printMatrix(double** matrix, int n, int m){
    
    cout << setw(8);
    cout << setprecision(5);
    cout << fixed;

    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
}