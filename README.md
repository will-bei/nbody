# nbody
Program that creates a 2000-frame animation of an n-body simulation

Description:
Program creates an n-body simulation. The dimensions of the simulation space is 4✕10^10 meters by 4✕10^10 meters. The gravitational constant used is 6.67✕10^10 m^3⋅kg^−1⋅s^−2. Initial masses of objects within the simulation are initialized randomly between 1✕10^22 kg to 1✕10^24 kg, inclusive.

Algorithms:
There are four main algorithms used in this program: one to sort the objects in order of decreasing mass, one to find and simulate collisions, one to calculate each object’s net acceleration, and one to find the position change of each object given a delta-time.

Sorting:
The purpose of sorting objects in order of decreasing mass is to simplify calculations in the next algorithm, which deals with collisions between objects. I used quicksort as my sorting algorithm.

Collisions:
The collisions algorithm is divided into two sections: the first is to check whether a collision has occurred, and the second is where calculations are made for the actual collision.

To check if any collisions have occurred, each object in the array checks if it is in contact with any object after it. If it has collided with another object, a perfectly inelastic collision is assumed. The resulting velocity is calculated using the formula (m1*v1 + m2*v2)/(m1+m2), where m1 and m2 are the first and second objects respectively and v1 and v2 are their respective velocities. The first object then increases in mass by the second object’s mass. The second object, which must be smaller than the first object due to the array being sorted in order of decreasing mass (see above), is moved to the end of the array, and the array size is decreased by one.

Note that each object stores and calculates its velocity by component form.

Acceleration:
To calculate each object’s net acceleration, individual accelerations by all other objects are added together. The formula used to calculate each acceleration is acceleration equals the gravitational constant times the second object’s mass divided by the distance between them squared. These calculations are repeated for all other objects in the array.

Note that each object stores and calculates its net acceleration by component form.

Position:
To calculate the each object’s change in position, each object’s velocity is changed by acceleration times the time step, or dt. Then, position is changed by its velocity times dt (because of the small value of dt, ½ a*dt2 is negligible).

Note that each object stores and calculates its position by component form.

Creating the .mp4 file:
To compile the program, enter g++ -o nbody_simulation nbody_3.3.cpp to command line.

The resulting executable program, nbody_simulation.exe, has three command line arguments: px, pz, and numberOfObjects.

px is the width of each frame in the mp4 file. pz is the height of each frame in the mp4 file. numberOfObjects is the number of MassObject objects the user wishes to have in the simulation.

As an example, to run the program with a frame size of 800 by 800 pixels and 300 objects, the following command would be entered: nbody_simulation.exe 800 800 300

(Note: on my own computer, runtime for these parameters takes approximately 25 minutes.)

After running the program, 2000 bmp files, each one a frame of the animation, will be created. Each will be named a five-digit integer (e.g. “00025.bmp”).

To create the actual animation, type the following into the command line:

magick mogrify -format jpeg *.bmp
ffmpeg -framerate 33 -i %05d.jpeg -s:v 600x400 -c:v libx264 -pix_fmt yuv420p animation.mp4

The result will be the animation file animation.mp4

An example result has been included in this repository.

ImageMagick must be installed for the animation to be created.
