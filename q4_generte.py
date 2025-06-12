# 生成 q4_patch.dat 文件内容（完整展开），共121节点，100个Q4单元
nx, ny = 10, 10  # 网格数量（单元数），节点数是 (nx+1)*(ny+1)
Lx, Ly = 1.0, 1.0
E = 1000
nu = 0.3
t = 1.0

dx = Lx / nx
dy = Ly / ny

nodes = []
node_id = 1
node_map = {}  # (i, j) -> node_id
for j in range(ny+1):
    for i in range(nx+1):
        x = round(i * dx, 5)
        y = round(j * dy, 5)
        # 边界条件：左边x固定，底部y固定
        bc_x = 1 if i == 0 else 0
        bc_y = 1 if j == 0 else 0
        nodes.append((node_id, bc_x, bc_y, x, y))
        node_map[(i, j)] = node_id
        node_id += 1

# Load: 节点在右上角 (nx, ny)，对x方向施加1N
load_node = node_map[(nx, ny)]
load_dir = 1  # x方向
load_val = 1.0

# 单元定义
elements = []
elem_id = 1
for j in range(ny):
    for i in range(nx):
        n1 = node_map[(i, j)]
        n2 = node_map[(i+1, j)]
        n3 = node_map[(i+1, j+1)]
        n4 = node_map[(i, j+1)]
        elements.append((elem_id, n1, n2, n3, n4, 1))
        elem_id += 1

# 生成文本内容
lines = []
lines.append("Q4_patch_10x10")
lines.append(f"{len(nodes)} 1 1 1")  # NUMNP, NUMEG, NLCASE, MODEX

# 节点数据
for nid, bx, by, x, y in nodes:
    lines.append(f"{nid} {bx} {by} {x:.5e} {y:.5e}")

# 荷载数据
lines.append("1 1")
lines.append(f"{load_node} {load_dir} {load_val:.5e}")

# 单元组数据
lines.append("2 1 1")  # ElementType=2(Q4), NUME, NUMMAT

# 单元数据
for eid, n1, n2, n3, n4, mat in elements:
    lines.append(f"{eid} {n1} {n2} {n3} {n4} {mat}")

# 材料数据
lines.append(f"1 {E:.5e} {nu:.5e} {t:.5e}")

from ace_tools import display_text_to_user
dat_content = "\n".join(lines)
display_text_to_user("q4_patch.dat", dat_content)

