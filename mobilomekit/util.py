import subprocess

def run_command(command, shell=False, timeout=None, input_text=None, env=None):
    """
    运行指定的命令并返回输出、错误和返回码。

    :param command: 要执行的命令（字符串或列表）。
    :param shell: 是否使用shell执行命令（默认为False）。
    :param timeout: 命令执行的超时时间（秒）。
    :param input_text: 传递给命令的标准输入（字符串）。
    :param env: 环境变量字典（默认为None，使用当前环境）。
    :return: 一个包含标准输出、标准错误和返回码的字典。
    """
    try:
        # 执行命令
        result = subprocess.run(
            command,
            shell=shell,
            timeout=timeout,
            input=input_text,
            capture_output=True,
            text=True,
            env=env
        )

        # 返回结果
        return {
            "stdout": result.stdout,
            "stderr": result.stderr,
            "returncode": result.returncode
        }
    except subprocess.TimeoutExpired:
        return {
            "stdout": "",
            "stderr": "Command timed out",
            "returncode": -1
        }
    except Exception as e:
        return {
            "stdout": "",
            "stderr": str(e),
            "returncode": -1
        }


# 示例用法
if __name__ == "__main__":
    # 调用系统命令
    result = run_command(["ls", "-l"])

    # 打印结果
    print("标准输出:", result["stdout"])
    print("标准错误:", result["stderr"])
    print("返回码:", result["returncode"])

    # 调用带输入的命令
    result = run_command(["grep", "hello"], input_text="hello world\nhi there\n")

    # 打印结果
    print("标准输出:", result["stdout"])
    print("标准错误:", result["stderr"])
    print("返回码:", result["returncode"])
    print(run_command.__annotations__)