{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Debugging"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Table of contents**<a id='toc0_'></a>    \n",
    "- 1. [Example](#toc1_)    \n",
    "- 2. [Numpy warnings](#toc2_)    \n",
    "- 3. [Scope bugs](#toc3_)    \n",
    "- 4. [Index bugs](#toc4_)    \n",
    "- 5. [Debugger](#toc5_)    \n",
    "  - 5.1. [In notebooks](#toc5_1_)    \n",
    "  - 5.2. [In modules](#toc5_2_)    \n",
    "\n",
    "<!-- vscode-jupyter-toc-config\n",
    "\tnumbering=true\n",
    "\tanchor=true\n",
    "\tflat=false\n",
    "\tminLevel=2\n",
    "\tmaxLevel=6\n",
    "\t/vscode-jupyter-toc-config -->\n",
    "<!-- THIS CELL WILL BE REPLACED ON TOC UPDATE. DO NOT WRITE YOUR TEXT IN THIS CELL -->"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 1,
   "metadata": {},
   "outputs": [],
   "source": [
    "import numpy as np"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Video:** [Why is a programming error called a bug?](https://www.youtube.com/watch?v=rhFSG-VyR_E)"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**General advice:**\n",
    "\n",
    "1. Code is always partly a **black box**\n",
    "2. **Print and plot results** to convince yourself (and others) that your results are sensible.\n",
    "3. Errors are typically something **very very simple**, look after that.\n",
    "4. If Python raises an error first try to **locate the line** where the error occurs.\n",
    "5. Your code can often run, but give you **unexpected behavior**."
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Most of the time spend programming is debugging!!**"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Assertions:** You can enforce assertations on e.g. numeric values"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 2,
   "metadata": {},
   "outputs": [],
   "source": [
    "x = -2\n",
    "y = x**2\n",
    "assert y > 0, f'x = {x}, y = {y}'"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Task:** Make the above assertion fail."
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Exceptions:** When code fails, it generates (*raises*) an exception "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 3,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "successful\n"
     ]
    }
   ],
   "source": [
    "x = 0.0\n",
    "try:\n",
    "    if x < 0.0:\n",
    "        raise ValueError('x cannot be negative')\n",
    "    else:\n",
    "        print('successful')\n",
    "except Exception as e:\n",
    "    print(e)"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Task:** Make it raise the exception and print *x cannot be negative*."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 1. <a id='toc1_'></a>[Example](#toc0_)\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Consider the following code:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 4,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "(3.340308817497993+0.5877852522924732j)\n"
     ]
    }
   ],
   "source": [
    "a = 0.8\n",
    "xlist = [-1,2,3]\n",
    "\n",
    "def myfun(xlist,a):\n",
    "    y = 0\n",
    "    for x in xlist:\n",
    "        z = x**a\n",
    "        y += z\n",
    "    return y\n",
    "\n",
    "y = myfun(xlist,a)\n",
    "print(y)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Problem:** Our result is a complex number. We did not expect that. Why does this problem arise?"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Find the error with print:**"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 5,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "x = -1 -> (-0.8090169943749473+0.5877852522924732j)\n",
      "x = 2 -> 1.7411011265922482\n",
      "x = 3 -> 2.4082246852806923\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "(3.340308817497993+0.5877852522924732j)"
      ]
     },
     "execution_count": 5,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "def myfun(xlist,a):\n",
    "    y = 0\n",
    "    for x in xlist:\n",
    "        z = x**a\n",
    "        print(f'x = {x} -> {z}') # temp\n",
    "        y += z\n",
    "    return y\n",
    "\n",
    "myfun(xlist,a)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Solution with an assert:**"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 6,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "z is not real for x = -1, but (-0.8090169943749473+0.5877852522924732j)\n"
     ]
    }
   ],
   "source": [
    "def myfun(xlist,a):\n",
    "    y = 0\n",
    "    for x in xlist:\n",
    "        z = x**a\n",
    "        assert np.isreal(z), f'z is not real for x = {x}, but {z}'\n",
    "        y += z\n",
    "    return y\n",
    "    \n",
    "try:\n",
    "    myfun(xlist,a)\n",
    "except Exception as e:\n",
    "    print(e)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Solution with if and raise exception:**"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 7,
   "metadata": {
    "tags": []
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "z is not real for x = -1, but (-0.8090169943749473+0.5877852522924732j)\n",
      "negative input number\n"
     ]
    }
   ],
   "source": [
    "def myfun(xlist,a):\n",
    "    y = 0\n",
    "    for x in xlist:\n",
    "        z = x**a\n",
    "        if not np.isreal(z):\n",
    "            print(f'z is not real for x = {x}, but {z}')\n",
    "            raise ValueError('negative input number') # an exception will be raised here  \n",
    "        y += z\n",
    "    return y\n",
    "\n",
    "try:\n",
    "    myfun(xlist,a)\n",
    "except Exception as e:\n",
    "    # we'll end up down here because the exception was raised. \n",
    "    print(e)"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Note:** You could also decide that the function should return e.g. **nan** when experiencing a complex number."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 8,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "nan"
      ]
     },
     "execution_count": 8,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "def myfun(xlist,a):\n",
    "    y = 0\n",
    "    for x in xlist:\n",
    "        z = x**a\n",
    "        if not np.isreal(z):\n",
    "            return np.nan\n",
    "        y += z\n",
    "    return y\n",
    "\n",
    "myfun(xlist,a)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 2. <a id='toc2_'></a>[Numpy warnings](#toc0_)"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Here we see an example of a *RuntimeWarning*."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 9,
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "C:\\Users\\gmf123\\AppData\\Local\\Temp\\ipykernel_11612\\3017935525.py:5: RuntimeWarning: invalid value encountered in log\n",
      "  y[i] = np.log(x)\n"
     ]
    },
    {
     "data": {
      "text/plain": [
       "array([       nan, 0.69314718, 1.09861229])"
      ]
     },
     "execution_count": 9,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "xlist = [-1,2,3]\n",
    "def f(xlist):\n",
    "    y = np.empty(len(xlist))\n",
    "    for i,x in enumerate(xlist):\n",
    "        y[i] = np.log(x)\n",
    "    return y\n",
    "\n",
    "f(xlist)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "You can **ignore all warnings**:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 10,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "array([       nan, 0.69314718, 1.09861229])"
      ]
     },
     "execution_count": 10,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "def f(xlist):\n",
    "    \n",
    "    y = np.empty(len(xlist))\n",
    "    for i,x in enumerate(xlist):\n",
    "        with np.errstate(all='ignore'):\n",
    "            y[i] = np.log(x)\n",
    "    return y\n",
    "\n",
    "f(xlist)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Better:** Decide what the code should do."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 11,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "array([      -inf, 0.69314718, 1.09861229])"
      ]
     },
     "execution_count": 11,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "def f(xlist):\n",
    "    \n",
    "    y = np.empty(len(xlist))\n",
    "    for i,x in enumerate(xlist):\n",
    "        if x <= 0:\n",
    "            y[i] = -np.inf\n",
    "        else:\n",
    "            y[i] = np.log(x)\n",
    "    return y\n",
    "\n",
    "f(xlist)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 3. <a id='toc3_'></a>[Scope bugs](#toc0_)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Global variables are dangerous:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 12,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "9.0\n"
     ]
    }
   ],
   "source": [
    "# a. define a function to multiply a variable with 5\n",
    "a = 5\n",
    "def f(x):\n",
    "    return a*x\n",
    "\n",
    "# many lines of code\n",
    "# many lines of code\n",
    "# many lines of code\n",
    "\n",
    "# z. setup the input and call f\n",
    "y = np.array([3,3])\n",
    "a = np.mean(y)\n",
    "b = np.mean(f(y))\n",
    "\n",
    "print(b)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Question:** What is the error?"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Conclusion:** \n",
    "\n",
    "1. Never use global variables, they can give poisonous **side effects**.\n",
    "2. Use a *positional argument* or a *keyword argument* instead. "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 4. <a id='toc4_'></a>[Index bugs](#toc0_)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 13,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "error found: index 10 is out of bounds for axis 0 with size 10\n"
     ]
    }
   ],
   "source": [
    "# a. setup\n",
    "N = 10\n",
    "x = np.linspace(1.3,8.2,N)\n",
    "y = 9.2\n",
    "\n",
    "# b. count all entries in x below y\n",
    "i = 0\n",
    "try:\n",
    "    while x[i] < y:\n",
    "        i += 1\n",
    "except Exception as e:\n",
    "    print(f'error found: {e}')"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Task:** Solve the problem."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## 5. <a id='toc5_'></a>[Debugger](#toc0_)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### 5.1. <a id='toc5_1_'></a>[In notebooks](#toc0_)"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Consider this example:"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "metadata": {},
   "outputs": [],
   "source": [
    "def f(x):\n",
    "    y = x-5\n",
    "    z = x\n",
    "    return np.log(y*z)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "nan\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "C:\\Users\\gmf123\\AppData\\Local\\Temp\\ipykernel_11612\\529491813.py:4: RuntimeWarning: invalid value encountered in log\n",
      "  return np.log(y*z)\n"
     ]
    }
   ],
   "source": [
    "x = 4\n",
    "q = f(x)\n",
    "print(q)"
   ]
  },
  {
   "attachments": {},
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "**Task:** Let us analyze why the problem arises with the *debugger*.\n",
    "\n",
    "1. Go to first line of `f(x)` two cells above. Press `F9` to create a *breakpoint*\n",
    "2. Go to the cell just above this cell. Press `Ctrl+Shift+Alt+Enter`\n",
    "3. Press `F10` to `Step Over`. Notice the value of y?\n",
    "4. Exit with `Shift+F5`\n",
    "\n",
    "**Extra:** Place the breakpoint at `x = 4`. Try out `F11` to `Step Into` the function `f(x)` and `Shift+F11` to `Step Out`\n",
    "\n",
    "**More details:** See [here](https://code.visualstudio.com/docs/datascience/jupyter-notebooks#_debug-cell)."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### 5.2. <a id='toc5_2_'></a>[In modules](#toc0_)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "We can also use the debugger in this case."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "metadata": {},
   "outputs": [],
   "source": [
    "%load_ext autoreload\n",
    "%autoreload 2"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "metadata": {},
   "outputs": [],
   "source": [
    "import mymodule"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 19,
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "nan\n"
     ]
    },
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "c:\\Users\\gmf123\\Dropbox\\Repositories\\course\\IntroProg-lectures\\1 - Debugging\\mymodule.py:6: RuntimeWarning: invalid value encountered in log\n",
      "  return np.log(y*z)\n"
     ]
    }
   ],
   "source": [
    "x = 4\n",
    "q = mymodule.g(x)\n",
    "print(q)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "base",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.15 (main, Nov 24 2022, 14:39:17) [MSC v.1916 64 bit (AMD64)]"
  },
  "toc-autonumbering": true,
  "vscode": {
   "interpreter": {
    "hash": "47ef90cdf3004d3f859f1fb202523c65c07ba7c22eefd261b181f4744e2d0403"
   }
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
