# `n=11` 单位圆盘覆盖证明包

本包给出定理

\[
r_{11}=0.379983853119838689722261606136093214308718956152215317\ldots
\]

的计算机辅助证明及独立重放程序。精确值定义见 `cover11_proof.md`：
`t_*` 是 `cover11_full_kkt_certificate.json` 所认证的 110 元有理二次
KKT 系统的唯一根的 `t` 坐标，`r_* = sqrt(t_*)`。

## 可信计算基

证明验证阶段只使用：

- Python 任意精度整数和 `fractions.Fraction`；
- Machin 公式的有理 `pi` 夹逼；
- 带显式 Lagrange 余项的定点 Taylor 正弦/余弦包围；
- 精确 Farkas 组合；
- 精确星—网/Kron 消元；
- 精确 `LDL^T` 正定性检查；
- 有理范数 Krawczyk 包含检查。

NumPy/SciPy/MPMath 只用于生成候选数据或被模块导入；所有决定性证书均由
有理验证路径重新检查。轮图权重另由
`cover11_wheel_row_validity_verify.py` 纯 `Fraction` 重放，不依赖生成时的
高精度浮点值。

## 环境

验证脚本目前使用绝对路径 `/mnt/data/<filename>`。最直接的运行方式是把
压缩包内容解压到 `/mnt/data`：

```bash
mkdir -p /mnt/data
tar -xzf cover11_complete_proof_bundle.tar.gz -C /mnt/data
cd /mnt/data
bash verify_all.sh quick
```

Python 依赖：Python 3.11+、NumPy、SciPy、MPMath。

## 验证模式

```bash
bash verify_all.sh quick
```

执行组合计数及 Brown 闭式计数核对、修正后的 Farkas 筛选、特殊残余轨道链接、
110 元 KKT 根、六组应力到 Kron 电导的精确重算、九维/八维局部坐标链接、
局部矩阵证书、候选证书结构、上界圆盘复形审计、轮图行有效性和全部分块结果聚合。

```bash
bash verify_all.sh full
```

在 `quick` 的基础上，从证书数据重新逐叶验证：

- 候选拓扑的 100,834 个直接叶盒；
- 5,302 个细分根盒中的 153,052 个外部叶盒和 1,953 个局部叶盒；
- 十边界轮图全部叶盒；
- 其余 53 个残余拓扑的全部 13,789 个叶盒。

## 主要文件

- `cover11_proof.md`：完整数学证明；
- `cover11_full_kkt_certificate.py/json`：精确候选根；
- `cover11_exact_upper_candidate_v2.py`：候选构型覆盖单位圆盘；
- `cover11_metric_Tplus.py/json`：3843 个组合类型的 Farkas 筛选；
- `cover11_topology_linkage_verify.py`：候选轨道 2547 与轮图轨道 1003 的精确链接；
- `cover11_candidate_*verify.py`：候选拓扑的全局和局部证书；
- `cover11_local_exact_audit2.py`：从六组 `weight_nums` 精确重算 Kron 电导；
- `cover11_residual_all_verify.py`：53 个非候选残余类型；
- `cover11_wheel_*`：十边界轮图；
- `cover11_complete_summary.json`：最终计数和最小裕量；
- `SHA256SUMS`：证明包文件哈希。

本证明尚未经过外部同行评审；包内结果表示所给数学归约和有限证书在当前
独立重放中全部通过。
