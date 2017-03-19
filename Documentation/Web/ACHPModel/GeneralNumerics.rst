********************
General Numerics
********************

One-dimensional Secant solver
=============================

If you have a function :math:`f`, which is a function of a single variable :math:`x` and the derivative of f cannot be found in a simple fashion, a secant solver can be used to find the value of :math:`x` that yields the equality :math:`f(x)=0`.  

The secant solver requires two initial guesses which set the search direction of the secant search.  It is not necessarily well-behaved, and can struggle to find a solution at times.  Also, if there are multiple solutions, you must start quite close to the solution that you want.  All that said, it is a simple model to implement, and performs quite admirably most of the time.  With the two initial guesses for :math:`x` of :math:`x_0` and :math:`x_1`, the function is evaluated at these points, yielding the functional values :math:`f_0=f(x_0)` and :math:`f_1=f(x_1)`.  The new guess for :math:`x` is then found from 

.. math::

    x_2=x_1-\dfrac{f_1}{\dfrac{f_1-f_0}{x_1-x_0}}
    
and the values of :math:`x_i` and :math:`f_i` are updated and this method is repeated until :math:`|f(x)|` is small enough.

The Scientific Python package includes an implementation of this method, and a simple example of its use is

.. ipython::

    #Import the newton function (implements secant method if no function derivative provided)
    In [1]: from scipy.optimize import newton
    
    # newton(function handle, x0), also see 
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.newton.html
    In [2]: newton(lambda x: x**2-0.4,1.5)
    
    # The exact solution - should be the same. 
    # a**b in Python is a to the power of b
    In [2]: 0.4**0.5

which provides a terribly inefficient solution for the square root of 0.4.

One-dimensional Bounded solver
==============================

If on the other hand, you know that the solution of :math:`f(x)=0` is bounded within some interval between a and b, a host of more powerful and robust solution methods are possible.

Bisection Method
----------------
The bisection method is guaranteed to converge so long as the range [a,b] brackets the solution and the function is continuous.  In the bisection method, the interval is repetitively chopped in half until the width of the interval :math:`w=|a-b|` is small enough.  First the values of the function are evaluated at :math:`a` and :math:`b`, and the midpoint of the interval :math:`m=(a+b)/2`, to yield the functional values of :math:`f(a)`, :math:`f(b)`, and :math:`f(m)`.  If :math:`f(a)f(m)>0`, :math:`a` and :math:`m` do not bracket the solution (solution must be between :math:`m` and :math:`b`), and the interval bounds are updated by

if :math:`f(a)f(m)>0`,

.. math::

    a=m
    
if :math:`f(a)f(m)<0`,

.. math::

    b=m
    
and the new midpoint is found from :math:`m=(a+b)/2`.  This method is applied until :math:`|a-b|` is small enough.  The same example as for the secant solver yields

.. ipython::

    #Import the bisect function (implements bisection method in 1-D)
    In [1]: from scipy.optimize import bisect
    
    # bisect(function handle, a,b), 
    # http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.bisect.html
    In [2]: bisect(lambda x: x**2-0.4,0,4)
    
    # The exact solution - should be the same. 
    # a**b in Python is a to the power of b
    In [2]: 0.4**0.5
    
Brent's Method
--------------
Brent's method, [#Brent]_ combines linear interpolation and quadratic interpolation with pure bisection if needed.  It is computationally efficient, robust, and stable, and used all over the code whenever a solution interval is known.  For well-behaved functions, the secant method is almost always faster, but for poorly behaved functions - like many in ACHP - Brent's method is required.

Python code given by::

    def Brent(f,a,b,macheps,t):
        """
        Using the algorithm from 
        
        Brent, R. P., Algorithms for Minimization Without Derivatives. 
        Englewood Cliffs, NJ: Prentice-Hall, 1973. Ch. 3-4.
        """
        fa=f(a)
        fb=f(b)
        c=a
        fc=fa
        if abs(fc)<abs(fb):
            # Goto ext: from Brent ALGOL code
            a=b
            b=c
            c=a
            fa=fb
            fb=fc
            fc=fa
        d=b-a
        e=b-a
        m=0.5*(c-b)
        tol=2*macheps*abs(b)+t
        while (abs(m)>tol and fb!=0):
            #See if a bisection is forced
            if abs(e)<tol or abs(fa) <= abs(fb):
                m=0.5*(c-b)
                d=e=m
            else:
                s=fb/fa
                if a==c:
                    #Linear interpolation
                    p=2*m*s
                    q=1-s
                else:
                    #Inverse quadratic interpolation
                    q=fa/fc
                    r=fb/fc
                    m=0.5*(c-b)
                    p=s*(2*m*q*(q-r)-(b-a)*(r-1))
                    q=(q-1)*(r-1)*(s-1)
                if p>0:
                    q=-q
                else:
                    p=-p
                s=e
                e=d
                m=0.5*(c-b)
                if 2*p<3*m*q-abs(tol*q) or p<abs(0.5*s*q):
                    d=p/q
                else:
                    m=0.5*(c-b)
                    d=e=m
            a=b
            fa=fb
            if abs(d)>tol:
                b+=d
            elif m>0:
                b+=tol
            else:
                b+=-tol
            fb=f(b)
            if fb*fc>0:
                # Goto int: from Brent ALGOL code
                c=a
                fc=fa
                d=e=b-a
            if abs(fc)<abs(fb):
                # Goto ext: from Brent ALGOL code
                a=b
                b=c
                c=a
                fa=fb
                fb=fc
                fc=fa
            m=0.5*(c-b)
            tol=2*macheps*abs(b)+t
        return b
    
.. _Numerical-Methods-NDsolve:

Multi-dimensional solver
========================
At the cycle-solver level, and elsewhere, is is common that multiple equations must be driven to zero by changing multiple parameters.  Beginning the consideration with the one-dimensional case, the Newton-Raphson method gives the solution for the next step from

.. math::

    x_{new}=x_{old}-\frac{f(x_{old})}{f'(x_{old})}
    
referring to the secant method above, you can see that the secant method is just the one-dimensional N-R method with a numerical approximation for the derivative.

In the multi-dimensional case, the problem to be solved is that to enforce the equality for the vector of nonlinear equations (within the convergence criterion), of

.. math::

    \mathbf f (\mathbf x)=0
    
where each of the functions :math:`f_1,f_2,...` can be functions of inputs :math:`x_1,x_2,...`.  The formulation of the problem in multidimensions is given by

.. math::

    \mathbf{x}_{new}=\mathbf{x}_{old}-\mathbf{J}'\mathbf{f}
    
where the Jacobian matrix :math:`\mathbf{J}` is given by

.. math::

    \mathbf{J}=\left[ \begin{array}{cccc} 
        \dfrac{\partial f_1}{\partial x_1} & \dfrac{\partial f_1}{\partial x_2} & \cdots& \dfrac{\partial f_1}{\partial x_n}\\
        \dfrac{\partial f_2}{\partial x_1} & \dfrac{\partial f_2}{\partial x_2} & \cdots& \dfrac{\partial f_2}{\partial x_n}\\
        \vdots & \vdots & \ddots & \vdots\\
        \dfrac{\partial f_n}{\partial x_1} & \dfrac{\partial f_n}{\partial x_2} & \cdots& \dfrac{\partial f_n}{\partial x_n}
        \end{array} \right]

and :math:`\mathbf{J}'` is the matrix inverse of :math:`\mathbf{J}`.  If the partial derivatives of the functions :math:`f_1,f_2,...` are known, they can be used directly in the calculation of the Jacobian matrix.  Otherwise, the Jacobian matix can be built with numeric derivatives.  The easiest way to build the numerical derivatives in the Jacobian matrix is to build the matrix by column.  If we call the vector of functional values at the iteration :math:`\hat {\mathbf{f}}`, each column is obtained by the formula 

.. math::

    \mathbf{J}_k=\frac{\partial \mathbf{f}}{\partial x_k}=\frac{\mathbf{f}(x_1,x_2,..x_k+\delta x,..., x_n)-\hat {\mathbf{f}}}{\delta x}
    
which forms the k-th column of :math:`\mathbf{J}`.  Only :math:`n+1` functional calls are required.  Scientific Python includes a slightly more advanced version of this algorithm with Jacobian updating, but the basic idea remains the same.  An example of this method is the set of equations 

.. math::

    x^2-2y-2=0
    
    x+y^2-1=0
    
which requires the python code

.. ipython::

    #Import the fsolve function (implements N-R multi-dimensional solve)
    In [1]: from scipy.optimize import fsolve
    
    # create an inline lambda function to return the values of the functions
    In [2]: func=lambda x: [x[0]**2-2*x[1]-2, x[0]+x[1]**2-1]
    
    #Actually solve the function    
    In [2]: x=fsolve(func,[0,0]); print(x)
    
    #Verify you have the right solution
    In [2]: func(x)

.. only:: html

    .. rubric:: References

.. [#Brent] Brent, R. P., Algorithms for Minimization Without Derivatives. Englewood Cliffs, NJ: Prentice-Hall, 1973. Ch. 3-4.  `Link to book <http://books.google.com/books?id=6Ay2biHG-GEC>`_