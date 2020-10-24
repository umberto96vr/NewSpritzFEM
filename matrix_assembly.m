syms Cg Ct Sg St real

Ry = [Ct 0 -St;
    0 1 0;
    St 0 Ct];
Rz = [Cg Sg 0;
    -Sg Cg 0;
    0 0 1];

T = Ry*Rz