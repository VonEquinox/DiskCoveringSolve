# 可选：证书发现脚本

本目录不进入证明接受链。接受只需根目录的 `python3 -S -B verify_all.py`。

这些脚本使用 NumPy、SciPy、mpmath；度量筛选还使用 numba。它们的成功标志和数值结果都不能替代精确验证。

先完成一次主重放，保留 `_runtime/` 的面表，然后运行：

```sh
python3 discovery/prepare_workspace.py
```

这会把发现脚本、相关核心模块及已有候选输入复制到 `_discovery_work/`，并链接/复制重建的面表。在该工作目录中可运行 `make_root15.py`、`make_local15.py`、`metric15.py`、`find_candidate15.py`、`build_candidate15.py`、`build_forest15.py`。数值发现可能因优化器版本、停止条件或浮点舍入而产生不同证书；只有精确验证器通过的证书才能被接受。

`force15.py` 用逐杆长度约束、锚点弧支持半平面及单位圆盘松弛作凸优化。取数值对偶力后，将其有理化，再通过生成树修正得到精确整数散度。所有候选排除均调用精确不等式检查，而不是信任求解器目标值。

本目录保留了最终候选的数值初值；不以完整保存所有随机启动、失败尝试或优化器日志作为可验证证明的前提。
