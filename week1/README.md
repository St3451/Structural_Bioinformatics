# Assignment 1

2PTC is a PDB ﬁle containing trypsin (chain E) and its inhibitor (chain I). The active site of trypsin includes residues 57, 102 and 195.

Download the PDB ﬁle and visualize the enzyme/inhibitor complex, showing the active site and the inhibitor residues (within 4 Å) around it. Write a fully automated pymol script to make the ﬁgure. Pymol scripts have the extension .pml and can be executed using "File -> Run Script". Make the picture as visually appealing, striking and clear as you can using pymol's functionality Hand in a single PDF ﬁle with  
* the image  
* a caption (description) for the image  
* the pymol script 


## Execution

###


### Caption

The image represent the trypsin serine protease and its inhibitor, they are contained in the 2PTC file from PDB.  
*	Pale yellow = trypsin enzyme  
*	Red = trypsin active site  
*	Green = trypsin inhibitor  
*	Blue = inhibitor residues 4 A° around the enzyme active site  
    *	Dark blue = atoms of the residues 4 A° around the active site  
    *	Light blue = atoms of the residues that are not 4 A° around the active site  					   
### Pymol script

delete all; fetch 2PTC  
set stick_ball; set stick_ball_ratio, 2  
bg_color gray; as ribbon  
select trp_enz, chain E; color paleyellow, trp_enz  
select trp_inhib, chain I; color forest, trp_inhib  
select act_site, trp_enz and (resi 57 or resi 102 or resi 195); color red, act_site; show sticks, act_site  
select inhib_close_aacids, (/2PTC/B/I/CYS\`14) or (/2PTC/B/I/LYS\`15) or (/2PTC/B/I/ALA`16); show sticks, inhib_close_aacids  
color skyblue, inhib_close_aacids  
select inhib_close_atoms, trp_inhib within 4 of act_site; color blue, inhib_close_atom;  
orient; zoom center, 30; set ray_opaque_background, on; ray 2400, 2000; save st3451.png  
