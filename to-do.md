# Notes from 5/1:

1. *Intro*: Motivation about how many results about quant gen and methods assume Gaussian distributions.
    Mention Levy walk macroevolution papers.
    Point out that usual quant gen has *no* sweeps, while stable distributions would have sweeps (since some s would be bigger than 1/N).
    Motivate with mixture model of effect sizes.
    If stable distributions are a better model, what does that tell us? We can better understand the *noise* in evolution, i.e., predictability and saltation-ness.
    (PETER)

2. Motivation: empirical SNPs distributions in real data (PETER)

3. Introducing infinitesimal model (in this case) (TODD)

4. Proof following BEV in stable case (TODD)

5. Theoretical results on stationary distributions (TODD)

6. Simulations (PETER)



# Notes from 4/4:

- Transformations: is the model well-defined?
    * this amounts to asking "if X, Y are indept draws from the trait distribution and xi is a a draw from seg variance, then if f = id
        the only f such that f((X+Y)/2 + xi) - (f(X) + f(Y))/2 is indept of f(X) and f(Y)?"
    * also note that GWAS effect sizes are estimated like E[f(X+\epsion) - f(X)].
    * but kinda clearly different underlying models give different behaviors for evolutionary trajectories and standing variation
        (e.g., are there sweeps?); so our question does make sense\
- Peter to look at the substitution process under adaptation and/or sweeps in genetic diversity, under a sudden change

# Notes from 3/21:

- add Cauchy mutations with Normal selection to simulations
- cancer codes look very nice in SNP plots - restrict to these?
- look at effect size distribution on 1+mu(x) since this is more analogous to cancer incidence?

# Notes from 3/14/23:

Peter to do simulation:
- mutation effect size Cauchy
- stabilizing selection with form in the writeup
- discretize time
- report distribution of offspring trait relative to parent mean (and trait distribution)

Todd to see if we can extend BEV reasonably.


# State of affairs, October 2019: We

- calculated decay exponents for SNP effects from Biobank: many are between 1 and 2.

- characterized fixed points in the absence of selection:

    * With Gaussian noise, fixed points are Cauchy (which is stable for the mean-ing process)
        plus Gaussian (the result of the noise).

    * With more general noise, we think that the fixed points are Cauchy (again, from the boundary)
        plus the distribution you get by doing $X_0 + (X_{10} + X_{11})/2 + (X_{20} + X_{21} + X_{22} + X_{23})/4 + \cdots$,
        where all $X$ are iid with the noise distribution.

- can explicitly solve for the Gaussian limit with Gaussian mortality selection and Gaussian noise;
    we suspect that Gaussian selection would destroy the Cauchy boundary contribution,
    so that this is the only solution.
    Note that a Cauchy reweighted by $x^{-k}$ has tails like $x^{(k+1)}$, so if $k \ge 1$ is in the Gaussian domain.

- think fecundity selection is more amenable to Fourier analysis than mortality selection,
    since under fecundity selection we multiply by a function that goes to zero.

- think that maybe we can say something about the tail decay given the tail decay of the noise and of selection.


To-do:

- Todd to see what he can say about solutions under stabilizing selection,
    and also directional selection: existence of Gaussian solutions? non-Gaussian solutions?

- Todd to also think more about fractal characteristic function:

    * conjecture: the set of all $2^k$-fractal characteristic functions as defined in the writeup
        are in correspondence with the "noise distributions" - given a distribution for $X$,
        this gives a fractal distribution by doing weighted tree averaging as above;
        and given such a distribution, the characteristic function of $X$ can be recovered
        by solving for it in the basic recursion.
        
- Peter or Todd to write down how fast adaptation goes in the general case by integrating the test function x
    against the mean measure evolution with fitness having positive slope.

- If it is obvious, saying how fast the measure drifts in neutral directions.

- Peter to make separate tail decay slope histograms for particular classes, eg. cancer-related codes

- Make up a mixture model that gets stable distributions of effect sizes.


Titles:

Abnormal Quantitative Genetics
Trait distributions under general noise processes
The infinitesimal model: beyond Normality.
Leaving Normal Behind.
A common framework for oligogenic and polygenic traits.


Tentative outline:

1. Motivation about how many results about quant gen and methods assume Gaussian distributions.
    Mention Levy walk macroevolution papers.
    Point out that usual quant gen has *no* sweeps, while stable distributions would have sweeps (since some s would be bigger than 1/N).
    Motivate with mixture model of effect sizes.
    If stable distributions are a better model, what does that tell us? We can better understand the *noise* in evolution, i.e., predictability and saltation-ness.

2. SNP effect distributions.

3. Introduction of infinitesimal model, and justification for it.

4. Solutions of stationary distribution without selection.
    Expressions for scale parameter under stable noise.

5. Speed of adaptation to directional selection.

6. Solutions of (or properties of them) to stabilizing selection,
    depending on form of the noise and form of selection.


Other notes from email:

Tentative intro: Quantitative genetics has a long and fairly successful history
of predicting what evolution, or at least artificial selection, is going to do,
and underpins lots of theoretical work on evolution, as well as phylogenetic
trait evolution. The assumption of Gaussian traits is almost universal, with
some notable exceptions; for instance, there's been some empirical support for
Levy-process models of trait evolution. The reasons for this are analytical
tractability and the CLT. Here, we look at non-Gaussian quantitative genetics,
showing that eg Cauchy trait distributions are stable limits, and that
empirical effect size distributions fall into the domain of attraction of
stable laws. This has the potential to reconcile empirical observations (that
much of the genetic variation in many traits is controlled by a few loci) and
quantitative genetics (that assumes implicitly that all loci have small
effect).

Todd:  I realized that there is a trivial (if contrived) example that gives
a Cauchy distribution for the stationary distribution under balancing
selection: if I assume Cauchy(0,a) noise and I take 1+mu(x) = (a+b)^2/a^2
* (b^2+x^2)/((a+b)^2+x^2), so that at phenotype x = 0, one has mortality one,
whereas mortality approaches a maximum of (a+b)^2/a^2 at infinity, then, unless
I’ve made a silly error (not unlikely) the population stationary distribution
is Cauchy(0,b).  (I’ve attached the relevant documents below so you can
remember what I’m talking about).

