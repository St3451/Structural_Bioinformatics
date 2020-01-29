delete all; fetch 2PTC
set stick_ball; set stick_ball_ratio, 2
bg_color gray; as ribbon
select trp_enz, chain E; color paleyellow, trp_enz
select trp_inhib, chain I; color forest, trp_inhib
select act_site, trp_enz and (resi 57 or resi 102 or resi 195); color red, act_site; show sticks, act_site
select inhib_close_aacids, (/2PTC/B/I/CYS`14) or (/2PTC/B/I/LYS`15) or (/2PTC/B/I/ALA`16); show sticks, inhib_close_aacids; color skyblue, inhib_close_aacids
select inhib_close_atoms, trp_inhib within 4 of act_site; color blue, inhib_close_atom;
orient; zoom center, 30; 
set ray_opaque_background, on; ray 2400, 2000
save stefano_pellegrini.png
