#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

#define UNROLLED_SIZE 9
#define INPUT_SIZE 3
#define KERNEL_SIZE 3
#define OUTPUT_SIZE 3

void printMatrix(char* name, int M, int N, double mat[M][N]) {
  printf("%s\n", name);

  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      printf ("%f ", mat[i][j]);
    }
    printf("\n");
  }
}

void flipMatrixHorizontal (int inputSize, double output[][inputSize], double input[][inputSize]) {
  for (int i =0; i< inputSize; i++) {
    for (int j=0; j< inputSize; j++) {
      output[inputSize - 1 - i][j] = input[i][j];
    }
  }
}

void flipMatrixVertical (int inputSize, double output[][inputSize], double input[][inputSize]) {
  for (int i =0; i< inputSize; i++) {
    for (int j=0; j< inputSize; j++) {
      output[i][inputSize - 1 - j] = input[i][j];
    }
  }
}

void unrollVectorToMatrix (int outputSize, double output[][outputSize], int inputSize, double input[inputSize]) {
  int count =0;
  for (int i =0; i<outputSize; i++) {
    for (int j=0; j<outputSize; j++) {
      output[i][j] = input[count];
      count++;
    }
  }
}

void unrollMatrixToVector(int outputSize, double unrolled[outputSize], int inputSize, double input[inputSize][inputSize]) {

  int count=0;
  for (int i=0; i<inputSize; i++) {
    for(int j=0; j<inputSize; j++) {
      unrolled[count] = input[i][j];
      count++;
    }
  }
}

void createCircularKernel(int outputSize, double unrolled[outputSize][outputSize], int inputSize, double input[inputSize]) {
  for (int i = 0; i < inputSize ; i++) {
    for (int j = 0; j < inputSize; j++) {
      int element = ((j-i)+inputSize)%inputSize;
      unrolled[j][i] = input[element];

    }
  }
}

void populateH (double H[9][9], double top[][3], double middle[][3], double bottom[][3], int offset) {
  for (int i=0;i<UNROLLED_SIZE;i++){
    for (int j=0; j<KERNEL_SIZE; j++) {
      if (i<KERNEL_SIZE) {
        H[i][j+offset]=top[i][j];
      }
      else if (i >= KERNEL_SIZE && i < KERNEL_SIZE*2) {
        H[i][j+offset] = middle[i-KERNEL_SIZE][j];
      }
      else {
        H[i][j+offset] = bottom[i-(KERNEL_SIZE*2)][j];
      }
    }
  }
}

void createMatrixh (int size, double h[][size], double kernel[][size], int rowNumber){
  double row[size];

  for (int i=0;i<size;i++) {
    row[i]= kernel[rowNumber][i];
  }

  for (int i = 0; i < size ; i++) {
    for (int j = 0; j < size; j++) {
      int element = ((j-i)+size)%size;
      h[j][i] = row[element];
    }
  }
}

int main()
{
  int i;
  int j;
  int m;
  int n;

  double input[3][3] = {{0, 1, 0},{1, 3, -1},{1, 2, 1}};
  double kernel[3][3] = {{0, 0, 0}, {1, 0 ,0 }, {1, -1, 0}};
  // 2d Convolution for Image Processing (Not CBLAS)
  double output[OUTPUT_SIZE][OUTPUT_SIZE];
  int sum = 0;
  for (m=0; m<OUTPUT_SIZE; m++) {
    for(n=0; n<OUTPUT_SIZE; n++) {
      int outerBound =  n - 1 + OUTPUT_SIZE;
      int innerBound = m - 1 + OUTPUT_SIZE;

      for (i= n-1 ; i<outerBound ; i++) {

        for (j = m-1; j<innerBound; j++) {

          if (i <0 || j < 0 || i > (OUTPUT_SIZE-1) || j > (OUTPUT_SIZE-1)) {
            sum += 0;
          }
          else {

            sum += input[j][i]*kernel[m-j+1][n-i+1];
          }
        }
      }
      output[m][n] = sum;
      sum = 0;
    }
  }
  printMatrix("Normal output", 3, 3, output);

  /*************************************************************
                        BLAS
  ***************************************************************/
  // Unroll Matrix (Done already above)
  // Make Kernel into circular toeplitz
  /************************1D*******************
  */
  /*
  double unrollKernel[UNROLLED_SIZE];
  unrollMatrixToVector(UNROLLED_SIZE, unrollKernel, KERNEL_SIZE, kernel);

  double circularKernel[UNROLLED_SIZE][UNROLLED_SIZE];
  createCircularKernel(UNROLLED_SIZE, circularKernel, UNROLLED_SIZE, unrollKernel);
*/

  /************************2D*******************************
  */
  // 2D convolution code for convolution of frequency

  // Construct h Matrix
  double h0[3][3];
  double h1[3][3];
  double h2[3][3];

  createMatrixh(KERNEL_SIZE, h2, kernel, 0);
  createMatrixh(KERNEL_SIZE,  h1, kernel, 1);
  createMatrixh (KERNEL_SIZE, h0, kernel, 2);

  // Construct H matrix
  double H[9][9];
  // Left most matricies
  populateH(H, h0, h1, h2, 0);
  // Middle matricies
  populateH(H, h2, h0, h1, KERNEL_SIZE);
  // Right most matricies
  populateH(H, h1, h2, h0, KERNEL_SIZE*2);

  // Create input matrix
  double unroll[UNROLLED_SIZE];
  double flipped[INPUT_SIZE][INPUT_SIZE];
  flipMatrixHorizontal(INPUT_SIZE, flipped, input);
  unrollMatrixToVector(UNROLLED_SIZE, unroll, INPUT_SIZE, flipped);

  // Do CBLAS stuff
  char no = 'N';
  double alpha = 1.0;
  double beta = 0.0;
  double y[9];

  cblas_dgemv(CblasRowMajor, CblasNoTrans, UNROLLED_SIZE, UNROLLED_SIZE, alpha, H, UNROLLED_SIZE, unroll, 1, beta, y, 1);

  double out[OUTPUT_SIZE][OUTPUT_SIZE];
  unrollVectorToMatrix(OUTPUT_SIZE, out, INPUT_SIZE, y);
  printMatrix("Cblas Output", 3, 3, out);

  return 0;
}