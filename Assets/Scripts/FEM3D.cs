using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using UnityEngine;
using System.Threading.Tasks;

namespace Assets.Scripts
{
    class FEM3D
    {
        public long E; // Young's modulus
        public float S; // Rod's section area

        public int[,] Links;
        public float[,] Coords_pre;
        public int[,] boundary_cond;
        public float[,] Concetrated_Forces;
        public int Num_of_timesteps;

        public float[,] U_total;
        public float[,] U_increment;
        public float[,] U_decrement;
        public float[,] Coords_current1;
        public float[,] Coords_current2;

        public int[] free_displ_coords;

        public FEM3D(int _Num_of_timesteps, int[,] _Links, float[,] _Coords_pre, int[,] _boundary_cond, float[,] _Concetrated_Forces)
        {
            E = 200000000000;
            S = 0.0001f;
            Num_of_timesteps = _Num_of_timesteps;
            Links = _Links;
            Coords_pre = _Coords_pre;
            boundary_cond = _boundary_cond;
            Concetrated_Forces = _Concetrated_Forces;
            
            Solver();

            U_increment = number_matrix_mult(U_total, (float)1 / Num_of_timesteps);
            U_decrement = number_matrix_mult(U_total, -(float)1 / Num_of_timesteps);

            Coords_current1 = Coords_pre;
        }

        // Finite Elements Method (FEM)
        void Solver()
        {
            // convertion of the concetrated forces
            float[,] Fc = new float[3 * (Coords_pre.GetUpperBound(1) + 1), 1];

            for (int i = 0; i < (Concetrated_Forces.GetUpperBound(1) + 1); i++)
            {
                int ind = (int)Concetrated_Forces[0, i];
                Fc[3 * ind, 0] = Concetrated_Forces[1, i]; //Fx
                Fc[3 * ind + 1, 0] = Concetrated_Forces[2, i]; //Fz
                Fc[3 * ind + 2, 0] = Concetrated_Forces[3, i]; //Fy
            }
            // creation of empty stiffness matrix
            int N = Coords_pre.GetUpperBound(1) + 1;
            float[,] K = new float[3 * N, 3 * N];

            // getting the coords of directions in nodes with non-zero displ 
            boundary_cond_converter(N);

            // stiffness matrix calculation
            for (int i = 0; i < Links.GetUpperBound(1) + 1; i++)
            {
                int[] ind = { Links[0, i], Links[1, i] };
                float[,] Ki = get_stiffness_matrix(get_float_matrix_part(Coords_pre, new int[] { 0, 1, 2 }, ind));
                Ki = transform_matrix(N, ind, Ki);
                K = matrix_sum(K, Ki);
            }

            // accounting for the boundary conditions
            float[,] K_new = get_float_matrix_part(K, free_displ_coords, free_displ_coords);
            float[,] F_new = get_float_matrix_part(Fc, free_displ_coords, new int[] { 0 });

            // solving the system of algebraic equations 
            float[,] U = new float[N * 3, 1];
            float[,] U_new = linsolveLU(K_new, F_new);
            int inh = 0;
            for (int i = 0; i < 3 * N; i++)
            {
                if (inh > free_displ_coords.GetUpperBound(0))
                {
                    U[i, 0] = 0;
                }
                else
                {
                    if (i == free_displ_coords[inh])
                    {
                        U[i, 0] = U_new[inh, 0];
                        inh = inh + 1;
                    }
                    else
                    {
                        U[i, 0] = 0;
                    }
                }
            }

            // getting the new coords of the nodes
            coords_converter(U);
        }

        // FEM auxiliary functions
        public void boundary_cond_converter(int N)
        {
            //
            int[] bc = new int[1];
            for (int i = 0; i < boundary_cond.GetUpperBound(1) + 1; i++)
            {
                int ind2 = boundary_cond[0, i];
                if (boundary_cond[1, i] == 1) // x
                {
                    int[] add = new int[1];
                    add[0] = 3 * ind2;
                    bc = bc.Concat(add).ToArray();
                }
                if (boundary_cond[2, i] == 1) // z
                {
                    int[] add = new int[1];
                    add[0] = 3 * ind2 + 1;
                    bc = bc.Concat(add).ToArray();
                }
                if (boundary_cond[3, i] == 1) // y
                {
                    int[] add = new int[1];
                    add[0] = 3 * ind2 + 2;
                    bc = bc.Concat(add).ToArray();
                }
            }
            bc = get_vector_part(bc, 1, bc.GetUpperBound(0));
            Array.Sort(bc, 0, bc.GetUpperBound(0) + 1);
            //
            free_displ_coords = new int[N * 3 - (bc.GetUpperBound(0) + 1)];

            int ind = 0;
            int ind_new = 0;
            for (int i = 0; i < N * 3; i++)
            {
                if (ind > bc.GetUpperBound(0))
                {
                    free_displ_coords[ind_new] = i;
                    ind_new = ind_new + 1;
                }
                else
                {
                    if (i == bc[ind])
                    {
                        ind = ind + 1;
                    }
                    else
                    {
                        free_displ_coords[ind_new] = i;
                        ind_new = ind_new + 1;
                    }
                }
            }
        }

        public float[,] get_stiffness_matrix(float[,] coords)
        {
            // Transform matrix definition
            float dX = coords[0, 1] - coords[0, 0]; // dX
            float dZ = coords[1, 1] - coords[1, 0]; // dZ
            float dY = coords[2, 1] - coords[2, 0]; // dY
            float[,] T = { { dX, dZ, dY, 0, 0, 0 }, { 0, 0, 0, dX, dZ, dY } };
            float L = (float)Math.Sqrt(dX * dX + dZ * dZ + dY * dY);
            T = number_matrix_mult(T, (float)1 / L);

            // ki definition
            float[,] ki = { { 1, -1 }, { -1, 1 } };
            ki = number_matrix_mult(ki, E * S / L);

            // Ki calculation
            float[,] Ki = matrix_matrix_mult(matrix_matrix_mult(matrix_transposition(T), ki), T);
            return Ki;
        }

        public float[,] transform_matrix(int N, int[] ind, float[,] Ki)
        {
            float[,] Ki_new = new float[N * 3, N * 3];

            Ki_new[3 * ind[0], 3 * ind[0]] = Ki[0, 0];
            Ki_new[3 * ind[0], 3 * ind[0] + 1] = Ki[0, 1];
            Ki_new[3 * ind[0], 3 * ind[0] + 2] = Ki[0, 2];
            Ki_new[3 * ind[0] + 1, 3 * ind[0]] = Ki[1, 0];
            Ki_new[3 * ind[0] + 1, 3 * ind[0] + 1] = Ki[1, 1];
            Ki_new[3 * ind[0] + 1, 3 * ind[0] + 2] = Ki[1, 2];
            Ki_new[3 * ind[0] + 2, 3 * ind[0]] = Ki[2, 0];
            Ki_new[3 * ind[0] + 2, 3 * ind[0] + 1] = Ki[2, 1];
            Ki_new[3 * ind[0] + 2, 3 * ind[0] + 2] = Ki[2, 2];

            Ki_new[3 * ind[1], 3 * ind[0]] = Ki[3, 0];
            Ki_new[3 * ind[1], 3 * ind[0] + 1] = Ki[3, 1];
            Ki_new[3 * ind[1], 3 * ind[0] + 2] = Ki[3, 2];
            Ki_new[3 * ind[1] + 1, 3 * ind[0]] = Ki[4, 0];
            Ki_new[3 * ind[1] + 1, 3 * ind[0] + 1] = Ki[4, 1];
            Ki_new[3 * ind[1] + 1, 3 * ind[0] + 2] = Ki[4, 2];
            Ki_new[3 * ind[1] + 2, 3 * ind[0]] = Ki[5, 0];
            Ki_new[3 * ind[1] + 2, 3 * ind[0] + 1] = Ki[5, 1];
            Ki_new[3 * ind[1] + 2, 3 * ind[0] + 2] = Ki[5, 2];

            Ki_new[3 * ind[0], 3 * ind[1]] = Ki[0, 3];
            Ki_new[3 * ind[0], 3 * ind[1] + 1] = Ki[0, 4];
            Ki_new[3 * ind[0], 3 * ind[1] + 2] = Ki[0, 5];
            Ki_new[3 * ind[0] + 1, 3 * ind[1]] = Ki[1, 3];
            Ki_new[3 * ind[0] + 1, 3 * ind[1] + 1] = Ki[1, 4];
            Ki_new[3 * ind[0] + 1, 3 * ind[1] + 2] = Ki[1, 5];
            Ki_new[3 * ind[0] + 2, 3 * ind[1]] = Ki[2, 3];
            Ki_new[3 * ind[0] + 2, 3 * ind[1] + 1] = Ki[2, 4];
            Ki_new[3 * ind[0] + 2, 3 * ind[1] + 2] = Ki[2, 5];

            Ki_new[3 * ind[1], 3 * ind[1]] = Ki[3, 3];
            Ki_new[3 * ind[1], 3 * ind[1] + 1] = Ki[3, 4];
            Ki_new[3 * ind[1], 3 * ind[1] + 2] = Ki[3, 5];
            Ki_new[3 * ind[1] + 1, 3 * ind[1]] = Ki[4, 3];
            Ki_new[3 * ind[1] + 1, 3 * ind[1] + 1] = Ki[4, 4];
            Ki_new[3 * ind[1] + 1, 3 * ind[1] + 2] = Ki[4, 5];
            Ki_new[3 * ind[1] + 2, 3 * ind[1]] = Ki[5, 3];
            Ki_new[3 * ind[1] + 2, 3 * ind[1] + 1] = Ki[5, 4];
            Ki_new[3 * ind[1] + 2, 3 * ind[1] + 2] = Ki[5, 5];

            return Ki_new;
        }

        public void coords_converter(float[,] U)
        {
            int N = (U.GetUpperBound(0) + 1) / 3;
            U_total = new float[3, N];
            int index = 0;

            for (int i = 0; i < N; i++)
            {
                U_total[0, index] = U[3 * i, 0];
                U_total[1, index] = U[3 * i + 1, 0];
                U_total[2, index] = U[3 * i + 2, 0];
                index = index + 1;
            }
        }


        // Operations with matrixes and vectors
        private float[,] get_float_matrix_part(float[,] A, int[] rows, int[] columns)
        {
            int NumRows = rows.GetUpperBound(0) + 1;
            int NumColumns = columns.GetUpperBound(0) + 1;
            float[,] B = new float[NumRows, NumColumns];
            for (int i = 0; i < NumRows; i++)
            {
                for (int j = 0; j < NumColumns; j++)
                {
                    B[i, j] = A[rows[i], columns[j]];
                }
            }
            return B;
        }

        private int[] get_vector_part(int[] A, int i1, int i2)
        {
            int[] B = new int[i2 - i1 + 1];
            for (int i = i1; i < i2 + 1; i++)
            {
                B[i - i1] = A[i];
            }
            return B;
        }

        static float[,] number_matrix_mult(float[,] A, float k)
        {
            float[,] B = new float[A.GetUpperBound(0) + 1, A.GetUpperBound(1) + 1];

            for (int i = 0; i < A.GetUpperBound(0) + 1; i++)
            {
                for (int j = 0; j < A.GetUpperBound(1) + 1; j++)
                {
                    B[i, j] = A[i, j] * k;
                }
            }
            return B;
        }

        private float[,] matrix_matrix_mult(float[,] A, float[,] B)
        {
            int Arows = A.GetUpperBound(0) + 1;
            int Brows = B.GetUpperBound(0) + 1;
            int Bcolumns = B.GetUpperBound(1) + 1;

            float[,] C = new float[Arows, Bcolumns];
            for (int i = 0; i < Arows; i++)
            {
                for (int j = 0; j < Bcolumns; j++)
                {
                    C[i, j] = 0;
                    for (int r = 0; r < Brows; r++)
                    {
                        C[i, j] += A[i, r] * B[r, j];
                    }
                }
            }
            return C;
        }

        private float[,] matrix_transposition(float[,] A)
        {
            int Arows = A.GetUpperBound(0) + 1;
            int Acolumns = A.GetUpperBound(1) + 1;
            float[,] AT = new float[Acolumns, Arows];
            for (int i = 0; i < Arows; i++)
            {
                for (int j = 0; j < Acolumns; j++)
                {
                    AT[j, i] = A[i, j];
                }
            }
            return AT;
        }

        private float[,] matrix_sum(float[,] A, float[,] B)
        {
            int rows = A.GetUpperBound(0) + 1;
            int columns = A.GetUpperBound(1) + 1;
            float[,] C = new float[rows, columns];
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    C[i, j] = A[i, j] + B[i, j];
                }
            }
            return C;
        }

        // Linsolve using LU decomposition
        private float[,] linsolveLU(float[,] A, float[,] b)
        {
            // LU Factorization
            int n = A.GetUpperBound(0) + 1;
            float[,] U = A;
            float[,] L = new float[n, n];


            for (int i = 0; i < n; i++)
                for (int j = i; j < n; j++)
                    L[j, i] = U[j, i] / U[i, i];

            for (int k = 1; k < n; k++)
            {
                for (int i = k - 1; i < n; i++)
                    for (int j = i; j < n; j++)
                        L[j, i] = U[j, i] / U[i, i];

                for (int i = k; i < n; i++)
                    for (int j = k - 1; j < n; j++)
                        U[i, j] = U[i, j] - L[i, k - 1] * U[k - 1, j];
            }

            // Ly = b -> y
            float[,] Y = new float[n, 1];
            Y[0, 0] = b[0, 0];
            for (int i = 1; i < n; i++)
            {
                Y[i, 0] = b[i, 0];
                for (int k = 0; k < i; k++)
                {
                    Y[i, 0] -= L[i, k] * Y[k, 0];
                }
            }

            // Ux = y -> x
            float[,] x = new float[n, 1];
            x[n - 1, 0] = Y[n - 1, 0] / U[n - 1, n - 1];
            for (int i = n - 2; i > -1; i--)
            {
                x[i, 0] = Y[i, 0];
                for (int k = n - 1; k > i; k--)
                {
                    x[i, 0] -= U[i, k] * x[k, 0];
                }
                x[i, 0] /= U[i, i];
            }

            return x;
        }


        // Update function
        public float[,] Update(int is_inverse)
        {
            if (is_inverse == 0)
            {
                Coords_current2 = matrix_sum(Coords_current1, U_increment);
            }
            else
            {
                Coords_current2 = matrix_sum(Coords_current1, U_decrement);
            }

            float[,] Rods_parameters = convert_parameters(Coords_current2);
            Coords_current1 = Coords_current2;
            return Rods_parameters;
        }

        // Function for converting the coordinates {Xstart, Ystart, Zstart, Xend, Yend, Zend} --> {Xc, Yc, Zc, Ay, Az, L}
        public float[,] convert_parameters(float[,] Coords)
        {
            float[,] Params = new float[7, Links.GetUpperBound(1) + 1];
            for (int i = 0; i < Links.GetUpperBound(1) + 1; i++)
            {
                int start = Links[0, i];
                int end = Links[1, i];

                Params[0, i] = (Coords[0, end] + Coords[0, start]) / 2; // centre coord x
                Params[1, i] = (Coords[1, end] + Coords[1, start]) / 2; // centre coord z 
                Params[2, i] = (Coords[2, end] + Coords[2, start]) / 2;  // centre coord y
                float dX = (float)Math.Abs(Coords[0, end] - Coords[0, start]);
                float dZ = (float)Math.Abs(Coords[1, end] - Coords[1, start]);
                float dY = (float)Math.Abs(Coords[2, end] - Coords[2, start]);
                float Ly = (float)Math.Sqrt(dZ * dZ + dX * dX);
                float L = (float)Math.Sqrt(dX * dX + dZ * dZ + dY * dY);
                Params[5, i] = L; // length

                Params[4, i] = Angle_calculation(Coords[0, start], Coords[1, start], Coords[0, end], Coords[1, end], Ly, dX, dZ, "Y", 0); // Ay
                Params[3, i] = Angle_calculation(Coords[0, start], Coords[2, start], Coords[0, end], Coords[2, end], L, Ly, dY, "Z", Coords[0, end] - Coords[0, start]); // Az
            }
            return Params;
        }

        public float Angle_calculation(float xs, float ys, float xe, float ye, float L, float dX, float dY, string Axis, float Normal)
        {
            float Angle;
            if ((dX < 0.0001) && (dY < 0.0001))
            {
                Angle = 0;
            }
            else
            {
                if (dX < 0.1)
                {
                    Angle = (float)Math.Acos(dX / L) / (float)Math.PI * 180;
                }
                else
                {
                    Angle = (float)Math.Asin(dY / L) / (float)Math.PI * 180;
                }
                if (xe > xs)
                {
                    if (ye > ys)
                    {
                        Angle = -Angle;
                    }
                }
                else
                {
                    if (ye > ys)
                    {
                        Angle = -180 + Angle;
                    }
                    else
                    {
                        Angle = 180 - Angle;
                    }
                }
            }

            if (Axis == "Z")
            {
                if (Normal > 0)
                {
                    Angle = -Angle;
                }
            }

            return Angle;
        }
    }
}
