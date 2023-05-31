import sys
import matplotlib.pyplot as plt
import scipy

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
    pwm = None # pwm[pos][A C G T] returns value at position of nucleotide
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
    pval = 0

    def __init__(self, name, alength, w, pwm):
        self.name = name
        self.w = w
        self.alength = alength
        self.pwm = pwm

    def graphThreshold(self):
        pwm_thresholds = []
        fig = plt.figure()
        fig.set_size_inches((10, 4))
        null_scores = self.rnd_scr_arr
        thresh = self.threshold
        ax = fig.add_subplot(1, 3,1)
        ax.hist(null_scores, bins=10)
        ax.axvline(x=thresh, color="red")
        ax.set_xlabel("Score")
        ax.set_ylabel("Frequency")
        ax.set_title(self.name)
        pwm_thresholds.append(thresh)
        fig.tight_layout()
        fig.show()
    
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
            temp_dist.pop(temp_dist.index(max(temp_dist)))
            i+=1
        thresh = max(temp_dist)
        self.threshold = thresh
        #self.graphThreshold()
        #print(str("SCORE ARRAY") + str(self.rnd_scr_arr))

    
    def setPval(self):
        peak_total = 0
        peak_motif = 0
        bg_total = 0

        bg_motif  = 0 #FIX THIS LATER
        print(len(self.rnd_scr_arr))
        peak_total = len(self.score_arr)
        for score in self.score_arr:
            if score >= self.threshold:
                peak_motif += 1
        for score in self.rnd_scr_arr:
            if score >= self.threshold:
                bg_motif += 1
        print("Peak total: " + str(peak_total) + " peak motif: " + str(peak_motif) + " Background total: " + str(peak_total) + " bgrnd motif: " + str(bg_motif))
        self.pval = KnownMotif.ComputeEnrichment(peak_total, peak_motif, peak_total, bg_motif)
    
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
        # list of scores. scores[i] should give the score of the substring sequence[i:i+n]. The default value for the score will be the maximum negative value so that the highest will always overwrite
        #scores = [0]*(len(sequence)-self.w + 1)
        scores = []
        i = 0
        # score the PWM on all spots in the sequence
        while i < ( len(sequence) - self.w + 1 ):
            scores.append(self.maxScore(sequence[i : i + self.w]))
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
            elif sequence[i] == 'N':
                revcomp = revcomp + 'N'
            i = i-1
        return revcomp

    def maxScore(self, sequence):
        """
        Score the sequence and its reverse complement and determine which is better
        """
        fwd_score = self.scoreSeq(sequence)
        bwd_score = self.scoreSeq(self.ReverseComplement(sequence))
        #print("Forward: " + str(fwd_score) + " Backward: " + str(bwd_score))
        if fwd_score >= bwd_score:
            return fwd_score
        else:
            return bwd_score

    def scoreSeq(self, sequence):
        """
        Given a sequence, find the score of that given sequence and return that score
        """
        score_temp = 0
        i = 0
        for nuc in str(sequence):
            # if the thing is an N, then the score is a 0
            if 'N' in nuc:
                score_temp += 0
            else:
                score_temp += float(self.pwm[self.nucs[nuc]][i])
            i += 1
        return score_temp
    

    
    def printName(self):
        print(self.name)

    def printPWM(self):
        nucs = ["A","C","G","T"]
        for row in range(len(self.pwm)):
            print(nucs[row], end = ": ")
            for letter in self.pwm[row]:
                print(str(letter)[0:4], end = " ")
            print()

    def printScoreArr(self):
        print(str(self.score_arr))
    
    def getThresh(self):
        return self.threshold

    def getScoreArr(self):
        return self.score_arr

    # To String method    
    def __str__(self):
        self.setPval() # MOVE THIS LINE
        return "MOTIF: " + self.name + " threshold score: " + str(self.threshold) + " pval:" + str(self.pval)