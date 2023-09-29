require(pwr)

# Whole sample
pwr.t.test(n = 470/2, 
           sig.level = 0.05, 
           power = 0.8, 
           type='two.sample', 
           alternative = 'greater')

# ketamine subset
pwr.t.test(n = 316/2, 
           sig.level = 0.05, 
           power = 0.8, 
           type='two.sample', 
           alternative = 'greater')
 
# esketamine subset
pwr.t.test(n = 82/2, 
           sig.level = 0.05, 
           power = 0.8, 
           type='two.sample', 
           alternative = 'greater')
