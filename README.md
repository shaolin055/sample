# sample 1
One of the fundamental tasks in my study was to find the local minimum positions for proton (H+) in an atomic structure, where a proton (H+) can be stable. Here I write a simple code to identify the initial possible locations for proton sites near oxygen (Protons typically bond to O with an O–H bond length of ∼1 Å in oxide materials). The workflow is as below:
- Identify symmetric oxygens in a structure. ( Once we can identify proton position around symmetric position, we can scale up to all oxygens)
- Mark 32 positions around each oxygen. ( 4 positions in each quadrant)
- If there is n oxygen in the atomic structure, we have n x 32 positions in the oxide structure.
- Prepare n x 32 input files for quantum computation (QC).
    We then used these files to relax the proton using QC and find where the proton ultimately be stable. Aggregating the last position gives us the possible proton sites in oxide materials. Identifying these sites helps us calculate activation energy for proton migration in oxide materials and identify which materials can easily incorporate proton.
