# Jesus Zeno Assignment 4

# Let's start by making dataframe of exon nucleotide frequency
exons = data.frame( first = c(18, 28, 33, 21),
                    internal = c(25, 28, 26, 21),
                    last = c(26, 22, 23, 29))
# Make dataframe of nucleotides
nucleotides_df = data.frame(nucleotide = c("A", "C", "G", "T"))
# Combine the nucleotides and exons df
exon_bp_df = cbind(nucleotides_df, exons) 
exon_bp_df

exon_matrix = data.matrix(exon_bp_df, as.character) #convert dataframe into a matrix
exon_matrix = exon_matrix[,-1] #remove the first column
exontable = as.table(exon_matrix) # convert matrix to a table

# Add the amount of specific nucleotides in all the exons. 
# Add the numbers in row up.
base_pairs_in_all_exons = margin.table(exontable,margin=1)
# Add the amount of base pairs in each exon. Add all numbers in column
base_pair_per_exon = margin.table(exontable,margin=2)
# Append these margins to the table.
summary_exon = addmargins(exontable)

# Let's make the barplot now
barplot(exontable, xlab="Exons", main="Nucleotides In Each Exon (ACGT)", 
        beside=TRUE)
# We can see that the barplot matches the DNA table we created earlier.
# In the first exon, we have G as the most nucleotides in the first exon, 
# C as the most in the internal exon, and T as the most in the last exon. 

#proportion of each base pair spread over the three exons
prop.table(exontable,margin=1) 
#proportion of base pairs in each exon
prop.table(exontable,margin=2)  

# Gives us the Chi-squared and p-values to see if there is any relation to 
# between the first exon and the internal exon. 
chisq.test(exon_bp_df$first,exon_bp_df$internal)
# P-value is high (>.05) so we reject the null hypothesis that the first exon
# and internal exon are independent. The results aren't significant.

# Gives us the Chi-squared and p-values to see if there is any relation to 
# between the first exon and the last exon. 
chisq.test(exon_bp_df$first,exon_bp_df$last)
# P-value is high (>.05) so we reject the null hypothesis that the first exon
# and last exon are independent. The results aren't significant.

# Gives us the Chi-squared and p-values to see if there is any relation to 
# between the first exon and the last exon. 
chisq.test(exon_bp_df$internal,exon_bp_df$last)
# P-value is high (>.05) so we reject the null hypothesis that the internal
# exon and last exon are independent. The results aren't significant.

# Conclusion: The p values are all high (>.05) so that tells me the different
# exons are not independent. We are rejecting the null hypothesis that
# the two variables in each instance are independent. This doesn't 
# necessarily mean that they are dependent on each other though. 


# Now we are going to do the exact same steps for the introns
introns = data.frame( first = c(28, 20, 22, 30),
                      internal = c(28, 20, 21, 31))
introns_bp_df = cbind(nucleotides_df, introns)
introns_bp_df
intron_matrix = data.matrix(introns_bp_df)
intron_matrix = intron_matrix[,-1]
introntable = as.table(intron_matrix)
introntable
bp_in_all_introns = margin.table(introntable,margin=1)
bp_per_intron = margin.table(introntable,margin=2)
summary_introns = addmargins(introntable)

barplot(introntable, xlab = "Introns", 
        main = "Nucleotides In Each Intron (ACGT)", beside = TRUE)
# This is interesting compared to the exon bar plot. The two introns seem
# to have almost identical amounts of each nucleotide as shown in the graph.  

prop.table(introntable,margin=1)
prop.table(introntable,margin=2)

# Gives us the Chi-squared and p-values to see if there is any relation to 
# between the two introns. 
chisq.test(introns_bp_df$first,introns_bp_df$internal)
# P-value is high (>.05) so we reject the null hypothesis that the first
# and internal introns are independent. The results aren't significant.

# Conclusion: The p values are all high (>.05) so that tells me the different
# introns are not independent. We are rejecting the null hypothesis that
# the two variables in each instance are independent. Again, this doesn't
# necessarily mean that they are dependent. 

# Part 2: Sampling

# Let's setup our infrastructure.
nucleotides = c('A','C','G','T') # Different Nucleotides
obs_first_exon = c(.18,.28,.33,.21) #first exon observed
obs_int_exon = c(.25,.28,.26,.21) #internal exon observed
obs_last_exon = c(.26,.22,.23,.29) #last exon observed
obs_first_intron = c(.28,.20,.22,.30) #first intron observed
obs_last_intron = c(.28,.20,.21,.31) #last intron observed

n=300 # Amount of nucleotides

# Initialize Counter and P-value for conformance
counter = 0
p=0

# Setup with alpha of .01  and p=0 and the while loop is <.05 
# since we are looking to see when it conforms.

# Look at first observed exon
while(p<.05){
    simulated_nucleotide = sample(nucleotides,n,replace=TRUE,prob=obs_first_exon)
    table(simulated_nucleotide)
    
    p = chisq.test(table(simulated_nucleotide), p=obs_first_exon)$p.value
    counter = counter + 1
  }
counter
table(simulated_nucleotide)
prop.table(table(simulated_nucleotide))
chisq.test(table(simulated_nucleotide), p=obs_first_exon)
# Chi-Square of 4.94 and p-value is >.05 (.175) so exit condition worked.

# Look at first internal exon
n=300 # Amount of nucleotides

# Initialize Counter and P-value for conformance
counter = 0
p=0
while(p<.05){
  simulated_nucleotide = sample(nucleotides,n,replace=TRUE,prob=obs_int_exon)
  table(simulated_nucleotide)
  
  p = chisq.test(table(simulated_nucleotide), p=obs_int_exon)$p.value
  counter = counter + 1
}
counter
table(simulated_nucleotide)
prop.table(table(simulated_nucleotide))
chisq.test(table(simulated_nucleotide), p=obs_int_exon)
# Chi-Square of .24 and p-value is >.05 (.969) so exit condition worked.

# Look at last exon observed
n=300 # Amount of nucleotides
# Initialize Counter and P-value for conformance
counter = 0
p=0
while(p<.05){
  simulated_nucleotide = sample(nucleotides,n,replace=TRUE,prob=obs_last_exon)
  table(simulated_nucleotide)
  
  p = chisq.test(table(simulated_nucleotide), p=obs_last_exon)$p.value
  counter = counter + 1
}
counter
table(simulated_nucleotide)
prop.table(table(simulated_nucleotide))
chisq.test(table(simulated_nucleotide), p=obs_last_exon)
# Chi-Square of 1.47 and p-value is >.05 (.687) so exit condition worked.

# Look at first intron observed
n=300 # Amount of nucleotides
# Initialize Counter and P-value for conformance
counter = 0
p=0
while(p<.05){
  simulated_nucleotide = sample(nucleotides,n,replace=TRUE,prob=obs_first_intron)
  table(simulated_nucleotide)
  
  p = chisq.test(table(simulated_nucleotide), p=obs_first_intron)$p.value
  counter = counter + 1
}
counter # First one to be 2 simulations instead of one.
table(simulated_nucleotide)
prop.table(table(simulated_nucleotide))
chisq.test(table(simulated_nucleotide), p=obs_first_intron)
# Chi-Square of 1.26 and p-value is >.05 (.737) so exit condition worked.

# Look at last intron observed
n=300 # Amount of nucleotides
# Initialize Counter and P-value for conformance
counter = 0
p=0
while(p<.05){
  simulated_nucleotide = sample(nucleotides,n,replace=TRUE,prob=obs_last_intron)
  table(simulated_nucleotide)
  
  p = chisq.test(table(simulated_nucleotide), p=obs_last_intron)$p.value
  counter = counter + 1
}
counter # Second one to be 2 simulations instead of one.
table(simulated_nucleotide)
prop.table(table(simulated_nucleotide))
chisq.test(table(simulated_nucleotide), p=obs_last_intron)
# Chi-Square of .138 and p-value is >.05 (.986) so exit condition worked. 


# Non-conforming part with alpha of .01
# For this we need to make the initial p-values 1 and then make the 
# while loop have p>.01. 


n=300 # Amount of nucleotides

# Initialize Counter and P-value for conformance
counter = 0
p=1

# Setup with alpha of .01  and p=0 and the while loop is <.05 
# since we are looking to see when it conforms.

# Look at first observed exon
while(p>.01){
  simulated_nucleotide = sample(nucleotides,n,replace=TRUE,prob=obs_first_exon)
  table(simulated_nucleotide)
  
  p = chisq.test(table(simulated_nucleotide), p=obs_first_exon)$p.value
  counter = counter + 1
}
counter # 135 simulations
table(simulated_nucleotide)
prop.table(table(simulated_nucleotide))
chisq.test(table(simulated_nucleotide), p=obs_first_exon)
# Chi-Square of 11.98 and p-value <.01 (.007) so the exit condition worked.


# Look at first internal exon
n=300 # Amount of nucleotides

# Initialize Counter and P-value for conformance
counter = 0
p=1
while(p>.01){
  simulated_nucleotide = sample(nucleotides,n,replace=TRUE,prob=obs_int_exon)
  table(simulated_nucleotide)
  
  p = chisq.test(table(simulated_nucleotide), p=obs_int_exon)$p.value
  counter = counter + 1
}
counter # 73 simulations
table(simulated_nucleotide)
prop.table(table(simulated_nucleotide))
chisq.test(table(simulated_nucleotide), p=obs_int_exon)
# Chi-Square of 12.83 and p-value <.01 (.005) so the exit condition worked.


# Look at last exon observed
n=300 # Amount of nucleotides
# Initialize Counter and P-value for conformance
counter = 0
p=1
while(p>.01){
  simulated_nucleotide = sample(nucleotides,n,replace=TRUE,prob=obs_last_exon)
  table(simulated_nucleotide)
  
  p = chisq.test(table(simulated_nucleotide), p=obs_last_exon)$p.value
  counter = counter + 1
}
counter # 70 simulations
table(simulated_nucleotide)
prop.table(table(simulated_nucleotide))
chisq.test(table(simulated_nucleotide), p=obs_last_exon)
# Chi-Square of 15.372 and p-value <.01 (.001) so the exit condition worked.


# Look at first intron observed
n=300 # Amount of nucleotides
# Initialize Counter and P-value for conformance
counter = 0
p=1
while(p>.01){
  simulated_nucleotide = sample(nucleotides,n,replace=TRUE,prob=obs_first_intron)
  table(simulated_nucleotide)
  
  p = chisq.test(table(simulated_nucleotide), p=obs_first_intron)$p.value
  counter = counter + 1
}
counter # 83 simulations
table(simulated_nucleotide)
prop.table(table(simulated_nucleotide))
chisq.test(table(simulated_nucleotide), p=obs_first_intron)
# Chi-Square of 14.229 and p-value <.01 (.002) so the exit condition worked.


# Look at last intron observed
n=300 # Amount of nucleotides
# Initialize Counter and P-value for conformance
counter = 0
p=1
while(p>.01){
  simulated_nucleotide = sample(nucleotides,n,replace=TRUE,prob=obs_last_intron)
  table(simulated_nucleotide)
  
  p = chisq.test(table(simulated_nucleotide), p=obs_last_intron)$p.value
  counter = counter + 1
}
counter # 45 simulations
table(simulated_nucleotide)
prop.table(table(simulated_nucleotide))
chisq.test(table(simulated_nucleotide), p=obs_last_intron)
# Chi-Square of 13.879 and p-value <.01 (.003) so the exit condition worked.