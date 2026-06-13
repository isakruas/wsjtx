TERN reference recordings (12 kHz mono, 16-bit PCM, signal at 1500 Hz,
realistic AWGN). Each file is one full T/R slot with the frame at a 1 s
lead-in.

  TERN_A_CQ_PU3IAR_GG40.wav   TERN Mode A, 30 s slot
  TERN_B_CQ_PU3IAR_GG40.wav   TERN Mode B, 60 s slot

Message: CQ PU3IAR GG40

Decode from the command line with the production decoder:
  jt9 --tern -b A -p 30 -f 1500 TERN_A_CQ_PU3IAR_GG40.wav
  jt9 --tern -b B -p 60 -f 1500 TERN_B_CQ_PU3IAR_GG40.wav
