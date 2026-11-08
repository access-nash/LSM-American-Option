# LSM-American-Option
I'm trying to value an American option using LSM algorithm as written in the original paper and comparing value using algorithm for the same in Yves Hilpisch book

A few main differences noted between the original LSM paper by Longstaff Schwartz and the algorithm in Yves Hilpisch book on Derivative Analytics are as follows-

i) The original paper uses only ITM paths whereas the book's algorithm uses all the paths. This is the major factor in the difference in valuations observed
ii) The basis function in the book uses just powers on the raw S, while the original paer uses weighted Laguerre polynomials. The latter is not something I was earlier familiar with, however, was easy to code from the original paper's algorithm
iii) The book version of the path simulation uses independent normals while the original LSM paper uses antithetic variates. 

The book's algorithm shows a valuation - "LSM American put value: 4.4041 (s.e. 0.0089)"
The LSM paper algorithm's implementation gives us - "LSM (paper) American put: 4.4822 (s.e. 0.0092)"

By using the same basis functions as in the book while only considering the in-the-money paths, I get a very similar valution as the original function - "LSM (paper) American put with Yves basis function: 4.4707 (s.e. 0.0092)". 

So we can conclude that the choice of basis functions or method used in path generation don't have a major impact on the valuation, however, considering only the ITM paths does lead to a higher estimate in the value by 0.07 
