using UnityEngine;
using System.Collections;

public class Class1 : MonoBehaviour
{
    public float Ax;
    public float Ay;
    public float Az;
    public float A2x;
    public float A2y;
    public float A2z;
    void Start()
    {
        transform.eulerAngles = new Vector3(0, 135,- 145);
        transform.position = new Vector3(0.5f, 0.5f, 0.5f);
        //transform.rotation = Quaternion.Euler(new Vector3(0, -145, 135));
        //transform.rotation = Quaternion.Euler(new Vector3(0, 215, 135));
    }
    void Update()
    {
        Ax = transform.eulerAngles.x;
        Az = transform.eulerAngles.z;
        Ay = transform.eulerAngles.y;       
        A2x = transform.localEulerAngles.x;
        A2z = transform.localEulerAngles.z;
        A2y = transform.localEulerAngles.y;  
    }
}