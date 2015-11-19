en = [0.9998 1.0002 0.9969 0.9955 1.0016 1.0018 1.0006 1.0016 0.9998 1.0009 ];
sn = [0.0821 0.1852 0.2404 0.2950 0.3544 0.4231 0.5084 0.6242 0.8113 1.4749 ];

eg = [0.9988 0.9987 0.9977 0.9984 1.0002 1.0005 1.0008 1.0013 1.0007 1.0000 ];
sg = [0.0821 0.1852 0.2404 0.2950 0.3544 0.4231 0.5084 0.6242 0.8113 1.4749 ];

plot(sn, en)
hold
plot(sg, eg)
legend('n', 'g')