using System.Collections;
using System.Collections.Generic;
using System.Text;
using System;
using UnityEngine;
using System.Globalization;

namespace Assets.Scripts
{
    public class SystemBehaviour : MonoBehaviour
    {
        // GameObjects
        public static GameObject[] Rods = new GameObject[42];
        public GameObject Rod;

        // System properties
        private float[,] Coords_pre;
        private int[,] Links;
        private float[,] Concetrated_Forces;
        private int[,] boundary_cond;

        // Number of timesteps
        public int Num_of_timesteps;

        // Physic model
        private FEM Model;

        // Timesteps
        public int timestep;
        public int is_inverse;

        // 
        public string Xcoords;
        public string Zcoords;
        public string Rods_start_nodes;
        public string Rods_end_nodes;
        public string BC_nodes;
        public string BC_X;
        public string BC_Z;
        public string Forces_nodes;
        public string Forces_X;
        public string Forces_Z;


        void Start()
        {
            var culture = (CultureInfo)CultureInfo.CurrentCulture.Clone();
            culture.NumberFormat.NumberDecimalSeparator = ".";

            // Coords_pre
            string[] Xcoords_vals = Xcoords.Split(' ');
            string[] Zcoords_vals = Zcoords.Split(' ');
            Coords_pre = new float[2, Xcoords_vals.GetUpperBound(0) + 1];
            for (int i = 0; i < Xcoords_vals.GetUpperBound(0) + 1; i++)
            {
                Coords_pre[0,i] = float.Parse(Xcoords_vals[i], culture);
                Coords_pre[1,i] = float.Parse(Zcoords_vals[i], culture);
            }

            // Links
            string[] Rods_start_nodes_vals = Rods_start_nodes.Split(' ');
            string[] Rods_end_nodes_vals = Rods_end_nodes.Split(' ');
            Links = new int[2, Rods_start_nodes_vals.GetUpperBound(0) + 1];
            for (int i = 0; i < Rods_start_nodes_vals.GetUpperBound(0) + 1; i++)
            {
                Links[0, i] = int.Parse(Rods_start_nodes_vals[i]);
                Links[1, i] = int.Parse(Rods_end_nodes_vals[i]);
            }

            // Boundary conditions
            string[] BC_nodes_vals = BC_nodes.Split(' ');
            string[] BC_X_vals = BC_X.Split(' ');
            string[] BC_Z_vals = BC_Z.Split(' ');
            boundary_cond = new int[3, BC_nodes_vals.GetUpperBound(0) + 1];
            for (int i = 0; i < BC_nodes_vals.GetUpperBound(0) + 1; i++)
            {
                boundary_cond[0, i] = int.Parse(BC_nodes_vals[i]);
                boundary_cond[1, i] = int.Parse(BC_X_vals[i]);
                boundary_cond[2, i] = int.Parse(BC_Z_vals[i]);
            }

            // Forces
            string[] Forces_nodes_vals = Forces_nodes.Split(' ');
            string[] Forces_X_vals = Forces_X.Split(' ');
            string[] Forces_Z_vals = Forces_Z.Split(' ');
            Concetrated_Forces = new float[3, Forces_nodes_vals.GetUpperBound(0) + 1];
            for (int i = 0; i < Forces_nodes_vals.GetUpperBound(0) + 1; i++)
            {
                Concetrated_Forces[0, i] = float.Parse(Forces_nodes_vals[i], culture);
                Concetrated_Forces[1, i] = float.Parse(Forces_X_vals[i], culture);
                Concetrated_Forces[2, i] = float.Parse(Forces_Z_vals[i], culture);
            }

            // FEM
            Model = new FEM(Num_of_timesteps, Links, Coords_pre, boundary_cond, Concetrated_Forces);

            // Objects coordinates: Xc, Zc, Angle_Y, L
            float[,] Params_pre = Model.convert_parameters(Coords_pre);
            for (int i = 0; i < Links.GetUpperBound(1) + 1; i++)
            {
                Rods[i] = Instantiate(Rod, new Vector3(Params_pre[0, i], 0, Params_pre[1, i]), Quaternion.Euler(new Vector3(0, Params_pre[2, i], 0)));
                Rods[i].transform.localScale = new Vector3(Params_pre[3, i], 0.025f, 0.025f); ;
            }

            // timestep initialization
            timestep = 0;
            is_inverse = 0;
        }

        void FixedUpdate()
        {
            if (is_inverse == 0)
            {
                timestep += 1;
            }
            else
            {
                timestep -= 1;
            }
            
            // Updating Object's coordinates and properties
            float[,] Rods_new_param = Model.Update(is_inverse);
            for (int i = 0; i < Links.GetUpperBound(1) + 1; i++)
            {
                Rods[i].transform.position = new Vector3(Rods_new_param[0, i], 0, Rods_new_param[1, i]);
                Rods[i].transform.rotation = Quaternion.Euler(new Vector3(0, Rods_new_param[2, i], 0));
                Rods[i].transform.localScale = new Vector3(Rods_new_param[3, i], 0.025f, 0.025f); ;
            }

            if (timestep == Num_of_timesteps)
            {
                is_inverse = 1;
            }
            if (timestep == 0)
            {
                is_inverse = 0;
            }
        }
    }

}
