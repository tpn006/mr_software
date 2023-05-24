import sys

class KnownMotif:
    forward_score = 0 # Running total of the highest score of each forward peak sequence
    score_arr = []
    pvalThresh = 0
    backward_score = 0 # Running total of the highest score of each reverse peak sequence
    random_score = 0 # Running total of the highest score of each random peak sequence
    rnd_scr_arr = []
    threshold = 0 # The score that divides matching from not matching
    name = "" # Name of the known motif
    w = 0 # Length of the known motif
    alength = 4 # number of nucleotides in the pwm
    pwm = [[0]*w]*alength # pwm[pos][A C G T] returns value at position of nucleotide
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
    forward_pval = 0
    backward_pval = 0

    def __init__(self, name, alength, w, pwm):
        self.name = name
        self.w = w
        self.alength = alength
        self.pwm = pwm
    def calcAndSetThreshold(self, background_sequences, pval):
        """
        Given a p-value threshold determine what score that is and set the threshold variable to that score
        """
        self.pvalThresh = pval
        for seq in background_sequences:
            self.scoreSequence(seq, "random")
        thresh = 0 # set this below to be the score threshold to obtain a p-value <0.01
        num = len(self.rnd_scr_arr) * pval
        temp_dist = self.rnd_scr_arr.copy()
        i = 1
        while i < num:
            print(temp_dist.pop(temp_dist.index(max(temp_dist))))
            i+=1
        thresh = max(temp_dist)
        self.threshold = thresh
        
    def ComputeEnrichment(peak_total, peak_motif, bg_total, bg_motif):
        """ 
        Compute fisher exact test to test whether motif enriched in bound sequences

        Parameters
        ----------
        peak_total : int
            Number of total peaks
        peak_motif : int
            Number of peaks matching the motif
        bg_total : int
            Number of background sequences
        bg_motif : int
            Number of background sequences matching the motif
       
        Returns
        -------
        pval : float
            Fisher Exact Test p-value    
        """
        pval = -1
        # your code here
        table = [[peak_motif, peak_total-peak_motif], [bg_motif, bg_total-bg_motif]]
        odds, pval = scipy.stats.fisher_exact(table)
        return pval

    def setPvalByFisher(self):
        #Assume threshold already set
        fwd_peak_pass = 0
        for seq in self.score_arr:
            if seq > self.threshold:
                fwd_peak_pass += 1  
        
        num_bg_pass= len(self.rnd_scr_arr) * self.pvalThresh

        self.pval = self.ComputeEnrichment(len(self.score_arr), fwd_peak_pass, len(self.rnd_scr_arr), num_bg_pass)

    def updateScore(self, val, which_score):
        """
        updates score by adding val to it
        """
        if "forward" in which_score:
            self.forward_score += float(val)
            self.score_arr.append(float(val))
        elif "random" in which_score:
            self.random_score += float(val)
            self.rnd_scr_arr.append(float(val))

    def scoreSequence(self, sequence, which_score):
        """
        Given a long sequence from a peak, find where this motif is best suited and then update the score.
        Calls a lot of scoreSeq with shorter sequences
        """
        scores = [0]*(len(sequence)-self.w + 1) # list of scores. scores[i] should give the score of the substring sequence[i:i+n]
        i = 0
        # score the PWM on all spots in the sequence
        while i < ( len(sequence) - self.w + 1 ):
            scores.append(self.scoreSeq(sequence[i : i + self.w]))
            i += 1
        scores.sort()
        if len(scores) != 0:
            self.updateScore(scores[len(scores)-1], which_score) #update score with the highest score found (last thing when sorted is the highest score)
    
    def ReverseComplement(self, sequence):
        """ Get the reverse complement of a sequence
                Parameters
                ----------
                sequence : str
                Sequence of nucleotides
                
                Returns
                -------
                revcomp : str
                Reverse complement of sequence
        """
        revcomp = ""
        i = len(sequence)-1
        while i >= 0:
            if sequence[i] == 'A':
                revcomp = revcomp + 'T'
            elif sequence[i] == 'C':
                revcomp = revcomp + 'G'
            elif sequence[i] == 'G':
                revcomp = revcomp + 'C'
            elif sequence[i] == 'T':
                revcomp = revcomp + 'A'
            i = i-1
        return revcomp

    def scoreSeq(self, sequence):
        """
        Given a sequence, find the score of that given sequence
        """
        score_fwd = 0
        i = 0
        for nuc in str(sequence):
            # if the thing is an N, then the score is a 0
            if 'N' in nuc:
                score_fwd += 0
            else:
                score_fwd += float(self.pwm[self.nucs[nuc]] [i])
            i += 1

        score_bwd = 0
        j = 0
        for nuc in self.ReverseComplement(str(sequence)):
            # if the thing is an N, then the score is a 0
            if 'N' in nuc:
                score_bwd += 0
            else:
                score_bwd += float(self.pwm[self.nucs[nuc]] [j])
                print(float(self.pwm[self.nucs[nuc]] [j]))
            j += 1

            
        if score_bwd < score_fwd:
            return score_fwd 
        else:
            return score_bwd 
    
    def printName(self):
        print(self.name)

    def printPWM(self):
        print(str(self.pwm))
    
    def getScore(self, which_score):
        if "forward" in which_score:
            return self.forward_score
        elif "backward" in which_score:
            return self.backward_score
        elif "random" in which_score:
            return self.random_score
    
    def printScore(self):
        print(str(self.score))

    def __str__(self):
        return "MOTIF: " + self.name + " threshold score: " + str(self.threshold) + " randoms score: " + str(self.rnd_scr_arr)