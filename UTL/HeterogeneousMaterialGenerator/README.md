It generates element wise heterogeneous material parameters for
   task 4, DECOVALEX 2019 by using the Gaussian distribution

run: hete_mat [mesh file name] [distribution file name]

Example of syntax of distribution file:
<data>
    <variable_name> porosity </variable_name>
    <gaussian> 0.13  0.01 0.11 0.14 </gaussian>
    <gaussian> 0.13  0.01 0.11 0.14 </gaussian>
    <gaussian> 0.13  0.01 0.11 0.14 </gaussian>
    <gaussian> 0.15  0.0276 0.097 0.185 </gaussian>
    <gaussian> 0.173  0.019 0.143 0.206 </gaussian>
    <gaussian> 0.193 0.029 0.15 0.249 </gaussian>
    <gaussian> 0.164  0.024 0.128 0.205 </gaussian>
    <gaussian> 0.1  0.001  0.095 0.11</gaussian>
</data>

The order of the data set of the distribution parameters:
  (mean, standard deviation, low range, high range)
is the same as the order of material IDs of mesh.

Note: there must be a space between the tag name and its values.
