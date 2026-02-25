echo "Spin-polarized current in Mn3Sn"
echo "============================================="
echo ""
echo "Command: symmetr res s.v E -f Mn3Sn.in"
echo ""
echo "==> Expected output:"
echo "
odd part of the response tensor:
X_ijk =
X_0pq =
⎡ 0    x₀₀₁  0⎤
⎢             ⎥
⎢x₀₁₀   0    0⎥
⎢             ⎥
⎣ 0     0    0⎦
X_1pq =
⎡x₁₀₀   0     0  ⎤
⎢                ⎥
⎢ 0    x₁₁₁   0  ⎥
⎢                ⎥
⎣ 0     0    x₁₂₂⎦
X_2pq =
⎡0   0     0  ⎤
⎢             ⎥
⎢0   0    x₂₁₂⎥
⎢             ⎥
⎣0  x₂₂₁   0  ⎦
even part of the response tensor:
X_ijk =
X_0pq =
⎡0   0     0  ⎤
⎢             ⎥
⎢0   0    x₀₁₂⎥
⎢             ⎥
⎣0  x₀₂₁   0  ⎦
X_1pq =
⎡ 0    0  x₁₀₂⎤
⎢             ⎥
⎢ 0    0   0  ⎥
⎢             ⎥
⎣x₁₂₀  0   0  ⎦
X_2pq =
⎡ 0    x₂₀₁  0⎤
⎢             ⎥
⎢x₂₁₀   0    0⎥
⎢             ⎥
⎣ 0     0    0⎦
"
echo "==> Output:"

symmetr res s.v E -f Mn3Sn.in

