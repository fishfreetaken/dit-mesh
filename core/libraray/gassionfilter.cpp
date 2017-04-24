typedef struct {  
    int cy;   
    int cx;  
} SIZE;  
  
typedef unsigned char *LPBYTE;  
  
typedef long LONG;  

//  一维高斯分布函数，用于平滑函数中生成的高斯滤波系数  
/* 
 *  @parameter sigma:     高斯函数参数 
 *  @parameter pdKernel:    高斯核函数模板 
 *  @parameter pnWidowSize:  高斯模板大小 
 */  
void CreatGauss(double sigma, double **pdKernel, int *pnWidowSize)  
{  
    LONG i;  
      
    int nCenter;//数组中心点  
    double dDis;//数组中一点到中心点距离  
  
    //中间变量  
    double dValue;  
    double dSum;  
    dSum = 0;  
      
    *pnWidowSize = 1+ 2*ceil(3*sigma);// [-3*sigma,3*sigma] 以内数据，会覆盖绝大部分滤波系数  
  
    nCenter = (*pnWidowSize)/2;  
    //生成高斯数据  
    for(i=0;i<(*pnWidowSize);i++)  
    {  
        dDis = (double)(i - nCenter);  
        dValue = exp(-(1/2)*dDis*dDis/(sigma*sigma))/(sqrt(2*3.1415926)*sigma);  
        (*pdKernel)[i] = dValue;  
        dSum+=dValue;  
    }  
    //归一化  
    for(i=0;i<(*pnWidowSize);i++)  
    {  
        (*pdKernel)[i]/=dSum;  
    }  
}  

//用高斯滤波器平滑原图像  
/* 
 *  @parameter sz   :  图像尺寸 
 *  @parameter pGray   :   图像灰度值 
 *  @parameter pResult:  图像 
 *  @parameter sigma:     高斯函数参数 
 */  
void GaussianSmooth(SIZE sz, LPBYTE pGray, LPBYTE pResult, double sigma)  
{  
    LONG x, y;  
    LONG i;  
      
    int nWindowSize;//高斯滤波器长度  
      
    int nLen;//窗口长度  
      
    double *pdKernel;//一维高斯滤波器  
      
    double dDotMul;//高斯系数与图像数据的点乘  
  
    double dWeightSum;//滤波系数总和  
  
    double *pdTemp;  
  
    nWindowSize = 1+ 2*ceil(3*sigma);// [-3*sigma,3*sigma] 以内数据，会覆盖绝大部分滤波系数  
  
  
    if ((pdTemp = (double *)malloc(sz.cx*sz.cy*sizeof(double)))==NULL)  
    {  
        printf("melloc memory for pdTemp failed!!");  
        exit(0);  
    }  
    if ((pdKernel = (double *)malloc(nWindowSize*sizeof(double)))==NULL)  
    {  
        printf("malloc memory for pdKernel,failed!!");  
        exit(0);  
    }  
  
    //产生一维高斯数据  
    CreatGauss(sigma, &pdKernel, &nWindowSize);  
  
    nLen = nWindowSize/2;  
  
    //x方向滤波  
    for(y=0;y<sz.cy;y++)  
    {  
        for(x=0;x<sz.cx;x++)  
        {  
            dDotMul = 0;  
            dWeightSum = 0;  
            for(i=(-nLen);i<=nLen;i++)  
            {  
                //判断是否在图像内部  
                if((i+x)>=0 && (i+x)<sz.cx)  
                {  
                    dDotMul+=(double)(pGray[y*sz.cx+(i+x)] * pdKernel[nLen+i]);  
                    dWeightSum += pdKernel[nLen+i];  
                }  
            }  
            pdTemp[y*sz.cx+x] = dDotMul/dWeightSum;     
        }  
    }  
  
    //y方向滤波  
    for(x=0; x<sz.cx;x++)  
    {  
        for(y=0; y<sz.cy; y++)  
        {  
            dDotMul = 0;  
            dWeightSum = 0;  
            for(i=(-nLen);i<=nLen;i++)  
            {  
                if((i+y)>=0 && (i+y)< sz.cy)  
                {  
                    dDotMul += (double)pdTemp[(y+i)*sz.cx+x]*pdKernel[nLen+i];  
                    dWeightSum += pdKernel[nLen+i];  
                }  
            }  
            pResult[y*sz.cx+x] = (unsigned char)dDotMul/dWeightSum;  
        }  
    }  
      
    free(pdTemp);//释放内存  
    free(pdKernel);  
}   