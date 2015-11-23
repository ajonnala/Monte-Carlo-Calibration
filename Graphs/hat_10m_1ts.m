% 1 time step
% 10 m particles
% hat kernel
% 50 strikes

eh = [0.3941 0.3942 0.3944 0.3945 0.3942 0.3933 0.3941 0.3926 0.3939 0.3939 0.3949 0.3933 0.3935 0.3947 0.3937 0.3945 0.3926 0.3943 0.3934 0.3942 0.3948 0.3942 0.3927 0.3935 0.3959 0.3939 0.3936 0.3938 0.3936 0.3945 0.3935 0.3940 0.3921 0.3940 0.3949 0.3938 0.3945 0.3949 0.3953 0.3941 0.3946 0.3938 0.3945 0.3940 0.3943 0.3941 0.3939 0.3933 0.3941 0.3943 ];
sh = [0.1232 0.2116 0.2341 0.2510 0.2653 0.2780 0.2896 0.3005 0.3107 0.3205 0.3300 0.3391 0.3481 0.3570 0.3657 0.3743 0.3829 0.3915 0.4001 0.4087 0.4174 0.4261 0.4350 0.4439 0.4530 0.4623 0.4718 0.4815 0.4914 0.5017 0.5122 0.5232 0.5345 0.5464 0.5588 0.5718 0.5855 0.6001 0.6156 0.6323 0.6504 0.6703 0.6923 0.7171 0.7454 0.7788 0.8194 0.8720 0.9474 1.0914 ];

eg = [0.9990 0.9980 0.9982 1.0009 1.0015 0.9998 0.9992 1.0002 1.0009 0.9997 0.9970 0.9954 0.9959 0.9972 0.9964 0.9941 0.9951 0.9978 0.9990 0.9985 0.9992 1.0031 1.0053 1.0024 1.0011 1.0028 1.0022 0.9996 0.9984 0.9998 0.9993 0.9982 0.9992 1.0021 1.0029 1.0010 1.0008 1.0022 1.0044 1.0046 1.0031 1.0011 0.9983 0.9967 0.9991 1.0007 1.0010 1.0020 1.0012 1.0003 ];
sg = [0.0492 0.1017 0.1209 0.1365 0.1502 0.1628 0.1748 0.1863 0.1974 0.2083 0.2190 0.2297 0.2403 0.2510 0.2617 0.2725 0.2834 0.2945 0.3057 0.3172 0.3289 0.3410 0.3533 0.3660 0.3791 0.3926 0.4067 0.4213 0.4365 0.4524 0.4690 0.4865 0.5050 0.5246 0.5455 0.5677 0.5916 0.6174 0.6454 0.6760 0.7100 0.7479 0.7909 0.8404 0.8988 0.9694 1.0588 1.1794 1.3627 1.7475 ];


plot(sg, eg, sh, eh)

% hold
% plot(stg2, eg2)
% hold
% plot(stg3, eg3)

legend('Gauss', 'Hat' )
title(' 1 Time step, NumPart = 10m, Strikes = 50')