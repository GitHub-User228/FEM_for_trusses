using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using UnityEngine;

namespace Assets.Scripts
{
    class FEM
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

        public FEM(int _Num_of_timesteps, int[,] _Links, float[,] _Coords_pre, int[,]  _boundary_cond, float[,] _Concetrated_Forces)
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
            float[,] Fc = new float[2 * (Coords_pre.GetUpperBound(1) + 1), 1];

            for (int i = 0; i < (Concetrated_Forces.GetUpperBound(1) + 1); i++)
            {
                int ind = (int)Concetrated_Forces[0, i];
                Fc[2 * ind, 0] = Concetrated_Forces[1, i];
                Fc[2 * ind + 1, 0] = Concetrated_Forces[2, i];
            }

            // creation of empty stiffness matrix
            int N = Coords_pre.GetUpperBound(1) + 1;
            float[,] K = new float[2 * N, 2 * N];

            // getting the coords of directions in nodes with non-zero displ 
            boundary_cond_converter(N);

            // stiffness matrix calculation
            for (int i = 0; i < Links.GetUpperBound(1) + 1; i++)
            {
                int[] ind = { Links[0, i], Links[1, i] };
                float[,] Ki = get_stiffness_matrix(get_float_matrix_part(Coords_pre, new int[] { 0, 1 }, ind));
                Ki = transform_matrix(N, ind, Ki);
                K = matrix_sum(K, Ki);
            }

            // accounting for the boundary conditions
            float[,] K_new = get_float_matrix_part(K, free_displ_coords, free_displ_coords);
            float[,] F_new = get_float_matrix_part(Fc, free_displ_coords, new int[] { 0 });

            // solving the system of algebraic equations 
            float[,] U = new float[N * 2, 1];
            float[,] U_new = linsolveLU(K_new, F_new);

            int inh = 0;
            for (int i = 0; i < 2 * N; i++)
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
                if (boundary_cond[1, i] == 1)
                {
                    int[] add = new int[1];
                    add[0] = 2 * ind2;
                    bc = bc.Concat(add).ToArray();
                }
                if (boundary_cond[2, i] == 1)
                {
                    int[] add = new int[1];
                    add[0] = 2 * ind2 + 1;
                    bc = bc.Concat(add).ToArray();
                }
            }
            bc = get_vector_part(bc, 1, bc.GetUpperBound(0));

            //
            free_displ_coords = new int[N * 2 - (bc.GetUpperBound(0) + 1)];

            int ind = 0;
            int ind_new = 0;
            for (int i = 0; i < N * 2; i++)
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
            float[,] T = new float[2, 4];
            float l_12 = coords[0, 1] - coords[0, 0];
            float m_12 = coords[1, 1] - coords[1, 0];
            float L = (float)Math.Sqrt(l_12 * l_12 + m_12 * m_12);
            T[0, 0] = l_12 / L;
            T[1, 2] = l_12 / L;
            T[0, 1] = m_12 / L;
            T[1, 3] = m_12 / L;

            // ki definition
            float[,] ki = { { 1, -1 }, { -1, 1 } };
            ki = number_matrix_mult(ki, E * S / L);

            // Ki calculation
            float[,] Ki = matrix_matrix_mult(matrix_matrix_mult(matrix_transposition(T), ki), T);

            return Ki;
        }

        public float[,] transform_matrix(int N, int[] ind, float[,] Ki)
        {
            float[,] Ki_new = new float[N * 2, N * 2];

            Ki_new[2 * ind[0], 2 * ind[0]] = Ki[0, 0];
            Ki_new[2 * ind[0], 2 * ind[0] + 1] = Ki[0, 1];
            Ki_new[2 * ind[0] + 1, 2 * ind[0]] = Ki[1, 0];
            Ki_new[2 * ind[0] + 1, 2 * ind[0] + 1] = Ki[1, 1];

            Ki_new[2 * ind[1], 2 * ind[0]] = Ki[2, 0];
            Ki_new[2 * ind[1], 2 * ind[0] + 1] = Ki[2, 1];
            Ki_new[2 * ind[1] + 1, 2 * ind[0]] = Ki[3, 0];
            Ki_new[2 * ind[1] + 1, 2 * ind[0] + 1] = Ki[3, 1];

            Ki_new[2 * ind[0], 2 * ind[1]] = Ki[0, 2];
            Ki_new[2 * ind[0], 2 * ind[1] + 1] = Ki[0, 3];
            Ki_new[2 * ind[0] + 1, 2 * ind[1]] = Ki[1, 2];
            Ki_new[2 * ind[0] + 1, 2 * ind[1] + 1] = Ki[1, 3];

            Ki_new[2 * ind[1], 2 * ind[1]] = Ki[2, 2];
            Ki_new[2 * ind[1], 2 * ind[1] + 1] = Ki[2, 3];
            Ki_new[2 * ind[1] + 1, 2 * ind[1]] = Ki[3, 2];
            Ki_new[2 * ind[1] + 1, 2 * ind[1] + 1] = Ki[3, 3];

            return Ki_new;
        }

        public void coords_converter(float[,] U)
        {
            int N = (U.GetUpperBound(0) + 1) / 2;
            U_total = new float[2, N];
            int index = 0;

            for (int i = 0; i < N; i++)
            {
                U_total[0, index] = U[2 * i, 0];
                U_total[1, index] = U[2 * i + 1, 0];
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

        // Function for converting the coordinates {Xstart, Ystart, Xend, Yend} --> {Xc, Yc, Angle_Y, L}
        public float[,] convert_parameters(float[,] Coords)
        {
            float[,] Params = new float[4, Links.GetUpperBound(1) + 1];
            for (int i = 0; i < Links.GetUpperBound(1) + 1; i++)
            {
                int start = Links[0, i];
                int end = Links[1, i];

                float[] A = angle(Coords[0, start], Coords[1, start], Coords[0, end], Coords[1, end]);

                Params[0, i] = (Coords[0, end] + Coords[0, start]) / 2; // centre coord x
                Params[1, i] = (Coords[1, end] + Coords[1, start]) / 2; // centre coord z 
                float dX = Math.Abs(Coords[0, end] - Coords[0, start]);
                float dY = Math.Abs(Coords[1, end] - Coords[1, start]);
                float L = (float)Math.Sqrt(Math.Pow(dX, 2) + Math.Pow(dY, 2));
                Params[3, i] = L; // length

                if (Coords[0, end] > Coords[0, start])
                {
                    if (Coords[1, end] > Coords[1, start])
                    {
                        if (dX < 0.1)
                        {
                            Params[2, i] = (-1) * (float)Math.Acos(dX / L) / (float)Math.PI * 180; // angle
                        }
                        else
                        {
                            Params[2, i] = (-1) * (float)Math.Asin(dY / L) / (float)Math.PI * 180; // angle
                        }
                    }
                    else
                    {
                        if (dX < 0.1)
                        {
                            Params[2, i] = (float)Math.Acos(dX / L) / (float)Math.PI * 180; // angle
                        }
                        else
                        {
                            Params[2, i] = (float)Math.Asin(dY / L) / (float)Math.PI * 180; // angle
                        }
                    }
                }
                else
                {
                    if (Coords[1, end] > Coords[1, start])
                    {
                        if (dX < 0.1)
                        {
                            Params[2, i] = (-1) * (90 + (float)Math.Asin(dX / L) / (float)Math.PI * 180); // angle
                        }
                        else
                        {
                            Params[2, i] = (-1) * (90 + (float)Math.Acos(dY / L) / (float)Math.PI * 180); // angle
                        }
                    }
                    else
                    {
                        if (dX < 0.1)
                        {
                            Params[2, i] = (90 + (float)Math.Asin(dX / L) / (float)Math.PI * 180); // angle
                        }
                        else
                        {
                            Params[2, i] = (90 + (float)Math.Acos(dY / L) / (float)Math.PI * 180); // angle
                        }
                    }
                }

            }
            return Params;
        }

        // Function to compute the Rod's Angle properly
        public float[] angle(float x1, float y1, float x2, float y2)
        {
            float[] A = new float[2];
            if (x2 > x1)
            {
                A[0] = 0;
                if (y2 > y1)
                {
                    A[1] = -1;
                }
                else
                {
                    A[1] = 1;
                }
            }
            else
            {
                A[0] = 90;
                if (y2 > y1)
                {
                    A[1] = -1;
                }
                else
                {
                    A[1] = 1;
                }
            }
            return A;
        }

    }

}
