Pressure Model Soft Body Simulation  
Tim Thirion  
CS 530 - Purdue University  
Fall 2004

The compile-time adjustable parameters are:

    static float    g_sphereRadius          - Adjust the initial radius of the sphere
    static uint32_t g_sphereSlices          - Adjust the number of slices of the sphere
    static uint32_t g_sphereStacks          - Adjust the number of stacks of the sphere
    static float    g_torusInnerRadius      - Adjust the inner radius of the torus
    static float    g_torusTubeRadius       - Adjust the tube radius of the torus
    static uint32_t g_torusThetaSlices      - Adjust the number of slices in the theta range [0, 2pi]
    static uint32_t g_torusPhiSlices        - Adjust the number of slices in the phi range [0, 2pi]
    static float    g_pointMass             - Adjust the mass of the point masses in the mesh
    static float    g_springConstant        - Adjust the spring constant for the sphere (the torus is different)
    static float    g_springDamping         - Adjust the damping constant for the springs
    static Vector3f g_gravity               - Adjust the strength and direction of gravity
    static float    g_internalPressure      - Adjust the initial internal pressure
    static float    g_targetPressure        - Adjust the target internal pressure (the mesh "inflates" to this value)
    static float    g_internalTemperature   - Adjust the internal temperature (and thus the buoyancy) of the mesh
    static Vector3f g_lightPosition         - Adjust the position of the light in the scene.

The rest of the relevant parameters should be adjustable at run-time. Read the menu for more information.

Enjoy!

-Tim

![Screenshot 1]
(https://raw.github.com/forscience/SoftBody/master/screenshots/screenshot1.png)
![Screenshot 2]
(https://raw.github.com/forscience/SoftBody/master/screenshots/screenshot2.png)
