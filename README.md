brune_python
============

Brune network synthesis in python

In order to synthesize a brune circuit, first put your impedance in the form

![](http://i.imgur.com/7MLnU0P.png)

i.e. a sum over poles with residues, and a constant offset. If this form is
not easily accessible, try fitting it numerically using
[vector fitting](http://github.com/PhilReinhold/vectfit_python).

Then use the function `vectfit_to_brune(poles, residues, offset)` to generate
a list of stages, and a final resistor. The stages are of the form
`(r1, l1, c2, l2, l3)`, coresponding to the brune model

![](http://i.imgur.com/idBRfby.png)

and the correspondence

![](http://i.imgur.com/El6CUxK.png)

Note that necessarily, either L1 or L3 will be negative. This would
seem unphysical except for equivalence to a set of coupled inductors.

Check that the model reproduces the data on
some test points with
`compose_brune_stages(1j*frequencies, stages, rfinal)`
