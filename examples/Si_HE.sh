echo "Hall effect in Si"
echo "============================================="
echo ""
echo "Command: symmetr res j E.B -f Si.in"
echo ""
echo "==> Expected output:"
echo "
odd part of the response tensor:
X_ijk =
X_0pq =
⎡0  0  0⎤
⎢       ⎥
⎢0  0  0⎥
⎢       ⎥
⎣0  0  0⎦
X_1pq =
⎡0  0  0⎤
⎢       ⎥
⎢0  0  0⎥
⎢       ⎥
⎣0  0  0⎦
X_2pq =
⎡0  0  0⎤
⎢       ⎥
⎢0  0  0⎥
⎢       ⎥
⎣0  0  0⎦
even part of the response tensor:
X_ijk =
X_0pq =
⎡0   0      0  ⎤
⎢              ⎥
⎢0   0    -x₂₁₀⎥
⎢              ⎥
⎣0  x₂₁₀    0  ⎦
X_1pq =
⎡  0    0  x₂₁₀⎤
⎢              ⎥
⎢  0    0   0  ⎥
⎢              ⎥
⎣-x₂₁₀  0   0  ⎦
X_2pq =
⎡ 0    -x₂₁₀  0⎤
⎢              ⎥
⎢x₂₁₀    0    0⎥
⎢              ⎥
⎣ 0      0    0⎦
"

echo "==> Output:"

symmetr res j E.B -f Si.in
