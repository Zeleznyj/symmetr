#!/bin/bash

echo "Conductivity and Anomalous Hall effect of Fe"
echo "============================================="
echo ""
echo "Command: symmetr res j E -f Fe.in"
echo ""
echo "==> Expected output:"
echo "
Applying Onsager relations.
even part of the response tensor:
⎡x₀₀   0    0 ⎤
⎢             ⎥
⎢ 0   x₀₀   0 ⎥
⎢             ⎥
⎣ 0    0   x₂₂⎦
odd part of the response tensor:
⎡ 0   -x₁₀  0⎤
⎢            ⎥
⎢x₁₀   0    0⎥
⎢            ⎥
⎣ 0    0    0⎦
"
echo "==> Output:"

symmetr res j E -f Fe.in


