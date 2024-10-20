import numpy as np
import scipy.stats as stats
import scipy.optimize as op

class Model:
    def __init__(self, parameters=None, coverage=100, readlength=100, times=None, labelingefficiency=0.05):
        """
        :param parameters: a list containing the transcript parameters [mu, nu, tau, lambda]
        :param coverage: Coverage, e.g. number of reads
        :param readlength: length of the reads, in nucleotides
        :param times: time series vector
        :param labelingefficiency: labeling efficiency. Needed for the generation of counts
        """
        # Initialize parameters
        if parameters is None:
            parameters = [1, 1, 1, 1]
        if times is None:
            times = [15, 30, 45, 60, 90, 120, 180]
        self.coverage = coverage  # coverage
        self.readlength = readlength  # read length
        self.times = np.array(times)  # time series
        self.labelingefficiency = labelingefficiency  # labeling efficiency
        self.mu, self.nu, self.tau, self.lam = parameters  # transcript parameters

        # Generate the new/total ratios given the transcript parameters
        self.nuc_total, self.nuc_unlabeled, self.nuc_labeled = self.nuclearfractions(
            mu=self.mu,
            nu=self.nu,
            tau=self.tau,
            t=self.times
        )
        self.cyt_total, self.cyt_unlabeled, self.cyt_labeled = self.cytosolicfractions(
            nuc_total=self.nuc_total,
            nu=self.nu,
            tau=self.tau,
            lam=self.lam,
            t=self.times
        )
        self.meas_nuc = self.nuclearratios(
            tau=self.tau,
            mu=self.mu,
            nu=self.nu,
            t=self.times
        ) # ratios for the nucleus
        self.meas_cyt = self.cytosolicratios(
            lam=self.lam,
            nuc_total=self.nuc_total,
            tau=self.tau,
            nu=self.nu,
            t=self.times
        )# ratios for the cytosol

        # Use the new/total ratios in order to simulate labeled and unlabeled counts.
        self.p_nuc = self.meas_nuc
        self.p_cyt = self.meas_cyt
        self.nuclear_labeledcounts = self.simulatecounts_nucleus()
        self.cytosolic_labeledcounts = self.simulatecounts_cytosol()
        self.nuc_newtotal_mlest, self.cyt_newtotal_mlest = self.perform_maximumlikelihood()

    def nuclearfractions(self, mu, nu, tau, t):
        """
        Computes the individual fractions of a transcript in the nucleus, namely the total amount, the unlabeled
        fraction, and the unlabeled fraction.

        :param mu: synthesis rate
        :param nu: nuclear degradation rate
        :param tau: nuclear export rate
        :param t: time series
        :return: total, unlabeled and labeled fraction in the nucleus
        """
        nuc_total = mu / (nu + tau)
        nuc_unlabeled = nuc_total * np.exp(-(nu + tau) * t)
        nuc_labeled = nuc_total - nuc_unlabeled

        return nuc_total, nuc_unlabeled, nuc_labeled

    def nuclearratios(self, tau, mu, nu, t):
        """
        Computes the nuclear new/total ratios based on the transcript parameters.
        Serves mainly as a wrapper for the nuclearfractions() function that is needed to specify an appropiate loss
        function below.

        :param mu: synthesis rate
        :param nu: nuclear degradation rate
        :param tau: nuclear export rate
        :param t: time series
        :return: new/total ratios in the nucleus at time point t of a time series
        """
        nuc_total, nuc_unlabeled, nuc_labeled = self.nuclearfractions(mu=mu, nu=nu, tau=tau, t=t)
        ratio = (nuc_labeled / nuc_total)

        return ratio

    def cytosolicfractions(self, nuc_total, nu, tau, lam, t):
        """
        Computes the individual fractions of a transcript in the cytosol, namely the total amount, the unlabeled
        fraction, and the unlabeled fraction.

        :param nuc_total: total fraction in the nucleus
        :param nu: nuclear degradation rate
        :param tau: nuclear export rate
        :param lam: cytosolic degradation rate
        :param t: time series
        :return:
        """
        cyt_total = nuc_total * tau / lam
        cyt_unlabeled = cyt_total * ((nu + tau) * np.exp(-lam * t) - lam * np.exp(-(nu + tau) * t)) / (nu + tau - lam)
        cyt_labeled = cyt_total - cyt_unlabeled

        return cyt_total, cyt_unlabeled, cyt_labeled

    def cytosolicratios(self, lam, nuc_total, nu, tau, t):
        """
        Computes the nuclear new/total ratios based on the transcript parameters.
        Serves mainly as a wrapper for the nuclearfractions() function that is needed to specify an appropiate loss
        function below.

        :param lam: cytosolic degradation rate
        :param nuc_total: total amount of RNA in the nucleus
        :param nu: nuclear degradation rate
        :param tau: nuclear export rate
        :param t: time series
        :return: new/total ratios in the cytosol at time point t of a time series
        """
        cyt_total, cyt_unlabeled, cyt_labeled = self.cytosolicfractions(
            nuc_total=nuc_total,
            nu=nu,
            tau=tau,
            lam=lam,
            t=t
        )
        ratio = (cyt_labeled / cyt_total)

        return ratio

    def simulatecounts_nucleus(self):
        """
        Simulate reads counts n based on the nuclear new/total ratios. The number of labeled counts n_lab are drawn
        from a binomial distribution Bin(N, (q * p(t))), where N is the coverage, p(t) the predicted new/total ratio
        at time point t and q the per-read labeling efficiency.

        :return: number of labeled reads n_lab and unlabeled reads n_unlab
        """
        n_lab = stats.binom.rvs(self.coverage, self.p_nuc)  # get number of labeled counts from binom distribution

        return n_lab

    def simulatecounts_cytosol(self):
        """
        Simulate reads counts n based on the cytosolic new/total ratios. The number of labeled counts n_lab are drawn
        from a binomial distribution Bin(N, (q * p(t))), where N is the coverage, p(t) the predicted new/total ratio
        at time point t and q the per-read labeling efficiency.

        :return: number of labeled reads n_lab and unlabeled reads n_unlab
        """
        n_lab = stats.binom.rvs(self.coverage, self.p_cyt)  # get number of labeled counts from binom distribution

        return n_lab


    def objective_globalefficiency(self, v, n_lab, coverage):
        """
        Specifies the objective function which is to be minimized/maximized using Maximum Likelihood.

        :param v: labeling efficiency which needs to be maximized by Maximum Likelihood
        :param n_lab: number of labeled counts
        :param coverage: coverage
        :param q: parameter that introduces variability in the estimation of the per-nucleotide labeling efficiency
        :return: probability of success to get n labeled counts from coverage N with labeling efficienc v * Q.
                 Result is returned as a negative number, since scipy uses a minimize function for Maximum Likelihood.
        """
        return - stats.binom.pmf(n_lab, coverage, v)

    def maximumlikelihood_newtotals(self, compartment):
        """
        Performs maximum likelihood for the nuclear labeled counts.

        :return: estimate for labeling efficieny for the nuclear labeled counts.
                 Equivalent to the 'real' new/total ratio
        """
        mlest = np.array([])
        for n in np.arange(len(compartment)):
            est = op.differential_evolution(
                func=self.objective_globalefficiency, # set predicted ratio as start,
                bounds=[[0, 1]],
                args=(compartment[n], self.coverage)
            )
            mlest = np.append(mlest, est.x[0])
        mlest = np.where(mlest < 1, mlest, 1)  # estimates < 0 become 0
        mlest = np.where(mlest > 0, mlest, 0)  # estimates > 1 become 1

        return mlest

    def perform_maximumlikelihood(self):
        """
        Performs maximum likelihood for the nuclear labeled counts.

        :return: estimate for labeling efficieny for the nuclear labeled counts.
                 Equivalent to the 'real' new/total ratio
        """
        mlest_nuc = self.maximumlikelihood_newtotals(self.nuclear_labeledcounts)
        mlest_cyt = self.maximumlikelihood_newtotals(self.cytosolic_labeledcounts)

        return mlest_nuc, mlest_cyt