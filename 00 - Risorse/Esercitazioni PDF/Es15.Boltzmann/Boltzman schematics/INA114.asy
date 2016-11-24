Version 3
SymbolType CELL
LINE Normal -5 -32 -5 -28
LINE Normal  5 -32  5 -28
LINE Normal -2 -25  2 -25
ARC Normal -5 -31 1 -25 -5 -28 -2 -25
ARC Normal -1 -31 5 -25  2 -25  5 -28
RECTANGLE Normal -32 -32 32 32
TEXT 0 0 CENTER 0 TI
WINDOW 0 0 -36 Center 0
WINDOW 3 0  36 Center 0
SYMATTR Prefix X
SYMATTR SpiceModel INA114.sub
SYMATTR Value INA114
SYMATTR Value2 INA114
SYMATTR Description Single Resistor Gain Programmable, Lowcost Instrumentation Amplifier

PIN -32 -24 LEFT 2
PINATTR PinName Rg1
PINATTR SpiceOrder 9

PIN -32  -8 LEFT 2
PINATTR PinName In-
PINATTR SpiceOrder 2

PIN -32   8 LEFT 2
PINATTR PinName In+
PINATTR SpiceOrder 1

PIN -32  24 LEFT 2
PINATTR PinName V-
PINATTR SpiceOrder 4

PIN  32  24 RIGHT 2
PINATTR PinName Ref
PINATTR SpiceOrder 8

PIN  32   8 RIGHT 2
PINATTR PinName Out
PINATTR SpiceOrder 5

PIN  32  -8 RIGHT 2
PINATTR PinName V+
PINATTR SpiceOrder 3

PIN  32 -24 RIGHT 2
PINATTR PinName Rg2
PINATTR SpiceOrder 10

PIN  0 32 BOTTOM 2
PINATTR PinName Fb
PINATTR SpiceOrder 6


