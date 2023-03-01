# Finite_element_method_application
Finite element method application used to calculate temperature distribution in an object

Whole program was written using C++. 

By using Finite element method (FEM) we can divide an object into smaller parts called elements, which form a mesh. By 
using a set of equations and approximating them we can, for example, calculate temperature distribution in an object.

The program reads all the important information about the object from .txt file and then calculates different matrixes and vectors for all of the elements of the object's mesh. After that it aggregates them and calculates temperature distribution. Results can then be put into a program such as ParaView to visualise flow of the heat. I used it to show the difference in temperature distribution in normal red brick (left) nd a brick used for building furnaces (right). 

![image](https://user-images.githubusercontent.com/78081338/222046153-5c0c7cbf-8a51-4fe3-9dde-a9eb9a347592.png)
