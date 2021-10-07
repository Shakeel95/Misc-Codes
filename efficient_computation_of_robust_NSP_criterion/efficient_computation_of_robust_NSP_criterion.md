## Efficient computation of robust NSP criterion

In the piecewise constant setting the (negative of the) robust NSP criterion can be shown to be shown to be weakly unimodal. Were it strictly unimodal, bracketing methods such as the golden section search could be used to find the minimum in O(\log n) time in the worst case. Unfortunately, since the criterion is only weakly unimodal constant regions can confuse the golden section search and we cannot guarantee to find the minimum.

A simple recursive extension to the golden section search (when a constant region is encountered) corrects this problem. The reclusive correction causes the algorithm to run in O(n) time in the worst case, however it still runs in O(\log n) time in the average case.

Both methods are implemented in R. Unfortunately recursive function calls are quite expensive in R, meaning the computational improvements are very modest in practice. 
