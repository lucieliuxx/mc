color Display Background white
display depthcue off
display projection Orthographic
axes location off
topo readvarxyz lamellarseq2.xyz
mol delete 0

mol addrep 1
mol modselect 1 1 name 1
mol modcolor 1 1 ColorID 23
mol modmaterial 1 1 Transparent
mol modstyle 1 1 VDW 0.4 12.0

mol addrep 1
mol modselect 2 1 name 2
mol modcolor 2 1 ColorID 1
mol modmaterial 2 1 AOChalky
mol modstyle 2 1 VDW 0.4 12.0

mol addrep 1
mol modselect 3 1 name 3
mol modcolor 3 1 ColorID 8
mol modmaterial 3 1 Transparent
mol modstyle 3 1 VDW 0.4 12.0

mol addrep 1
mol modselect 4 1 name 4
mol modcolor 4 1 ColorID 15
mol modmaterial 4 1 AOChalky
mol modstyle 4 1 VDW 0.4 12.0

mol addrep 1
mol modselect 5 1 name 5
mol modcolor 5 1 ColorID 32
mol modmaterial 5 1 AOChalky
mol modstyle 5 1 VDW 0.4 12.0
mol delrep 0 1
