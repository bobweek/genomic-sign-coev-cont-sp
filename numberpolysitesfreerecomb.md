
# Number of polymorphic sites in a population with freely recombining loci

Since each gene tree is identically distributed, they have the same expected time to most recent common ancestor; $\mathbb ET_{MRCA}\approx 4N_e$. With a per-lineage mutation rate of $\mu$ (so that the expected number of mutations in the population per generation is $2N_e\mu$), the expected total number of mutations since the most recent common ancestor (not only the ones that make it into a sample or land on an extant lineage), is then $(4N_e\times2N_e)\mu=8N_e^2\mu$. This provides the upper bound for the expected number of polymorphic sites in a sample. The expected time interval between the arrival of each mutation is $1/2N_e\mu$. Labelling the mutations by the order they appear going backwards in time, the expected time for the $\ell$th mutation is $\ell/2N_e\mu$. The probability that the $\ell$th mutation corresponds to a polymorphic sites is equal to the probability that it lands on an extant lineage. The number of extant lineages, $j$, $t$ units of time ago is the maximum $j$ that satisfies

$$2N_e\sum_{i=j+1}^{2N_e}2/i(i-1)<t.$$

Write this maximal $j$ as $j(t)$. Hence, the probability that the $\ell$th mutation corresponds to a polymorphic site is $1/j(\ell/2N_e\mu)$. Note that

$$2N_e\sum_{i=j+1}^{2N_e}2/i(i-1)=4N_e\left(\frac{2N_e-1}{2N_e}-\frac{j-1}{j}\right)\approx \frac{4N_e}{j}.$$

So $j(\ell/2N_e\mu)$ approximately satisfies

$$\frac{1}{j}<\frac{\ell}{8N_e^2\mu}\implies j(\ell/2N_e\mu)\approx \mathrm{floor}\left(\frac{8N_e^2\mu}{\ell}\right).$$

Randomly choosing mutations, mutation $\ell$ is chosen with probability $1/8N_e^2\mu$. Then, out of all the sites where a mutation has landed since the most recent common ancestor, the probability that a randomly chosen site is currently polymorphic is

$$p=\frac{1}{8N_e^2\mu}\sum_{\ell=1}^{8N_e^2\mu}\frac{1}{j(\ell/2N_e\mu)}\approx\frac{1}{64N_e^4\mu^2}\sum_{\ell=1}^{8N_e^2\mu}\ell=\frac{8N_e^2\mu+1}{16N_e^2\mu}.$$

When $N_e\mu\gg1$, $p\approx1/2$ (which seems weird). Then, my approximation of the expected number of polymorphic sites in a population is $S_B=4N_e^2\mu$. To see how similar this is to Watterson's result $S_W=4N_e\mu(\gamma+\ln(2N_e))$  ($\gamma$ is Euler's constant) in the case of no recombination, we can numerically evaluate. Plugging in $N_e=1e4$ and $\mu=1e-4$, I get $S_B\approx4e4$ and $S_W\approx42$, vastly different quantities...
