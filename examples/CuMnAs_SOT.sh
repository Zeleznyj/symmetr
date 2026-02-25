echo "Projected spin-orbit torque in CuMnAs"
echo "============================================="
echo ""
echo "Command: symmetr res s E -f CuMnAs.in -p 5 -p2 6"
echo ""
echo "==> Expected output:"
echo "
even part of the response tensor:
⎡ 0   x₀₁  0⎤
⎢           ⎥
⎢x₁₀   0   0⎥
⎢           ⎥
⎣ 0    0   0⎦
odd part of the response tensor:
⎡ 0   0  x₀₂⎤
⎢           ⎥
⎢ 0   0   0 ⎥
⎢           ⎥
⎣x₂₀  0   0 ⎦
First part of the response tensor, atom 6
⎡ 0    -x₀₁  0⎤
⎢             ⎥
⎢-x₁₀   0    0⎥
⎢             ⎥
⎣ 0     0    0⎦
Second part of the response tensor, atom 6
⎡ 0   0  x₀₂⎤
⎢           ⎥
⎢ 0   0   0 ⎥
⎢           ⎥
⎣x₂₀  0   0 ⎦
"
echo "==> Output:"

symmetr res s E -f CuMnAs.in -p 5 -p2 6
