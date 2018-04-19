function[C1,C2]=Mute(P1,P2,k)
A1=floor(P1/(2^k));
B1=P1-A1*(2^k);
A2=floor(P2/(2^k));
B2=P2-A2*(2^k);
C1=A1*(2^k)+B2;
C2=A2*(2^k)+B1;

