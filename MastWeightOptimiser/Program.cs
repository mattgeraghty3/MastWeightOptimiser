using System;
using System.Threading;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using RES.CSVReader;
using System.IO;
using Clerc2011;
using System.Globalization;
using Microsoft.SolverFoundation.Common;
using Microsoft.SolverFoundation.Services;

namespace MastWeightOptimizer
{
    class Program
    {
        private static readonly CultureInfo _culture = new CultureInfo("en-GB");

        static void Main(string[] args)

        {
/****** To Change the Stack Size of the .exe change the value of the stackSize Variable***********************************************************************************************************************************************************/
            const int stackSize = 419430400;
            Thread th = new Thread(() =>
            {
/**************************************************************************************************************************************************************************************************************************************************/
/************************Part 1: Declare and Initialize all Variables and Read in Values from CSVs*************************************************************************************************************************************************/
                CSVReader mastReader = new CSVReader(new FileStream("MastLocations.csv", FileMode.Open));

                var value = int.Parse(mastReader["MastID"][0]);

                var masts = new Mast[mastReader.RowCount];
                var mastIDs = new int[mastReader.RowCount];

                for (int i = 0; i < mastReader.RowCount; i++)
                {
                    masts[i] = new Mast(double.Parse(mastReader["MastX"][i], _culture), double.Parse(mastReader["MastY"][i], _culture));
                    mastIDs[i] = int.Parse(mastReader["MastID"][i]);
                }

                CSVReader turbineReader = new CSVReader(new FileStream("TurbineLocationsAndMastWeights.csv", FileMode.Open));

                var turbines = new Turbine[turbineReader.RowCount];
                var turbineIDs = new int[turbineReader.RowCount];

                var Weight = new double[turbineReader.RowCount, mastReader.RowCount];

                for (int i = 0; i < turbineReader.RowCount; i++)

                {
                    turbines[i] = new Turbine(double.Parse(turbineReader["TurbineX"][i], _culture), double.Parse(turbineReader["TurbineY"][i], _culture));
                    turbineIDs[i] = int.Parse(turbineReader["TurbineID"][i]);

                    for (int j = 0; j < mastReader.RowCount; j++)
                    {
                        Weight[i, j] = double.Parse(turbineReader["Mast" + (j + 1) + "Weight"][i]);
                    }

                }

                var _binTopoUncertaintyCalculator = new BinTopoUncertaintyCalculator(new CoefficientOfVariationCalculator(new MastToTurbineVariationCalculator(0.1, 1000), new SpeedupVariationCalculator(0.5)), new CorrelationCalculator(new AngleCorrelationCalculator(), new DistanceCorrelationCalculator(1000)));

                var binReader = new CSVReader(new FileStream("BinInformation.csv", FileMode.Open));

                int count = binReader.RowCount;
                int number = turbineReader.RowCount * mastReader.RowCount;

                double[] variance = new double[count];

                Term total = 0;
                Term Energy = 0;
                Term Uncertainty = 0;

                int[] T = new int[binReader.RowCount];
                int[] M = new int[binReader.RowCount];

/**********************************Part 2: Create Solver Model, Decisions and Constraints **********************************************************************************/

                SolverContext context = SolverContext.GetContext();
                Model model = context.CreateModel();

                BinInformation[] bins = new BinInformation[binReader.RowCount];

                Decision[,] MWeight = new Decision[turbineReader.RowCount, mastReader.RowCount];

                for (int row = 0; row < turbineReader.RowCount; row++)
                {
                    for (int col = 0; col < mastReader.RowCount; col++)

                    {
                        MWeight[row, col] = new Decision(Domain.RealNonnegative, "MastWeight" + (row+1) + (col+1));
                        model.AddDecision(MWeight[row, col]);
                        model.AddConstraint("limit" + row + col, 0 <= MWeight[row, col] <= 1);
                        model.AddConstraint("nonzero" + row + col, MWeight[row, col] >= 0);
                        /*** Remove this line if no initial values desired *********************/
                        MWeight[row, col].SetInitialValue(Weight[row, col]);
                        /***********************************************************************/
                    }
                }

                /******* Sum of mast Weights must equal 1 constraint****************************/
                for (int i = 0; i < turbineReader.RowCount; i++)
                {
                    Term balance = 0;
                    Term holder = 0;

                    for (int j = 0; j < mastReader.RowCount; j++)
                    {
                        holder = MWeight[i, j];
                        balance += holder;
                    }
                    model.AddConstraint("Balances" + i, balance == 1);
                }
                /********************************************************************************/

/************************************Part 3: Read in Bin Information from CSV**********************************************************************************************************/

                for (int i = 0; i < binReader.RowCount; i++)
                {
                    int mastID = int.Parse(binReader["MastID"][i]);
                    int turbineID = int.Parse(binReader["TurbineID"][i]);

                    int mastIndex = GetIndex(mastIDs, mastID);
                    int turbineIndex = GetIndex(turbineIDs, turbineID);
                    bins[i].BinEnergy = double.Parse(binReader["Energy"][i], _culture);
                    bins[i].Direction = double.Parse(binReader["Direction"][i], _culture);
                    bins[i].SpeedUp = double.Parse(binReader["SpeedUp"][i], _culture);
                    bins[i].SensitivityFactor = double.Parse(binReader["Sensitivity"][i], _culture);
                    bins[i].Turbine = turbines[turbineIndex];
                    bins[i].Mast = masts[mastIndex];
                    T[i] = int.Parse(binReader["TurbineID"][i], _culture);
                    M[i] = int.Parse(binReader["MastID"][i], _culture);
                }

/***************************Part 4: Define the Objective Function: I.e. Uncertainty************************************************************************************************************/

                double[,] MastWeights = new double[turbineReader.RowCount, mastReader.RowCount];
                for (int i = 0; i < count; i++)
                {
                    int mast1 = M[i] - 1;
                    int turbine1 = T[i] - 1;

                    if (bins[i].BinEnergy > 0 && bins[i].SensitivityFactor > 0)
                    {
                        variance[i] = _binTopoUncertaintyCalculator.CalculateVariance(bins[i]);

                        if (variance[i] > 0)
                        {
                            for (int j = 0; j <= i; j++)
                            {
                                int mast2 = M[j] - 1;
                                int turbine2 = T[j] - 1;

                                if (bins[j].BinEnergy > 0 && bins[j].SensitivityFactor > 0 && variance[j] > 0)
                                {

                                    var singleResult = _binTopoUncertaintyCalculator.Calculate(bins[i], bins[j], variance[i], variance[j]) * bins[i].BinEnergy * bins[j].BinEnergy * ((MWeight[turbine1, mast1])) * ((MWeight[turbine2, mast2]));

                                    if (i != j)
                                        singleResult *= 2;

                                    total += singleResult;
                                }
                            }
                        }
                    }
                    Energy += bins[i].BinEnergy * ((MWeight[turbine1, mast1]));

                }

                total = Model.Sqrt(total);
              
                Uncertainty = total / Energy;
/**********************************Part 5: Define the solving algorithm and its parameters*********************************************************************************************************/

                HybridLocalSearchDirective MastWeightOptimizer = new HybridLocalSearchDirective();

                /*******Update these solver parameters by changing these variables******************/
                MastWeightOptimizer.RunUntilTimeout = false;
                //MastWeightOptimizer.TimeLimit = 60000*60*10;
                MastWeightOptimizer.EqualityTolerance = 0.001;
                MastWeightOptimizer.WaitLimit = -1;
                /*************************************************************************************/

/**************************************Part 6: Solve the Model ********************************************************************************************************************/

                model.AddGoal("Uncertainty", GoalKind.Minimize, Uncertainty);

                Solution solution = context.Solve(MastWeightOptimizer);

                Report report = solution.GetReport();

/**********************************Part 7: Print Report to the Console *******************************************************************************************************/

                Console.WriteLine("{0}", report);

/********************Part 8: Convert solution to double floating points and check calculation**********************************************************************************/

                for (int i = 0; i < turbineReader.RowCount; i++)
                {
                    for (int j = 0; j < mastReader.RowCount; j++)
                    {
                        MastWeights[i, j] = ((MWeight[i, j].ToDouble()));   
                    }
                }

                double FinalEnergy = 0;
                double FinalSigma = 0;
                double FinalUncertainty = 0;

                var DubBinTopoUncertaintyCalculator = new BinTopoUncertaintyCalculator(new CoefficientOfVariationCalculator(new MastToTurbineVariationCalculator(0.1, 1000), new SpeedupVariationCalculator(0.5)), new CorrelationCalculator(new AngleCorrelationCalculator(), new DistanceCorrelationCalculator(1000)));

                for (int i = 0; i < count; i++)
                {
                    int mast1 = M[i] - 1;
                    int turbine1 = T[i] - 1;

                    if (bins[i].BinEnergy > 0 && bins[i].SensitivityFactor > 0)
                    {
                        variance[i] = DubBinTopoUncertaintyCalculator.CalculateVariance(bins[i]);

                        if (variance[i] > 0)
                        {
                            for (int j = 0; j <= i; j++)
                            {
                                int mast2 = M[j] - 1;
                                int turbine2 = T[j] - 1;

                                if (bins[j].BinEnergy > 0 && bins[j].SensitivityFactor > 0 && variance[j] > 0)
                                {

                                    double singleResult = DubBinTopoUncertaintyCalculator.Calculate(bins[i], bins[j], variance[i], variance[j]) * bins[i].BinEnergy * bins[j].BinEnergy * MastWeights[turbine1, mast1] * MastWeights[turbine2, mast2];

                                    if (i != j)
                                        singleResult *= 2;

                                    FinalSigma += singleResult;
                                }
                            }
                        }
                    }

                    FinalEnergy += bins[i].BinEnergy * MastWeights[turbine1, mast1];
                }

                FinalSigma = Math.Sqrt(FinalSigma);

                FinalUncertainty = FinalSigma / FinalEnergy;

/************************Part 9: Print the final results for Uncertainty, Energy, Sigma and Mast Weights to the Console***********************************************************************************************/

                Console.WriteLine("************* Optimised Mast-To-Turbine Weihghtings**********************************************");
                Console.WriteLine("");
                Console.Write("Turbine#     ");
                for (int i = 0; i<mastReader.RowCount; i++)
                {
                    Console.Write("MastWeight{0}    ", (i + 1));
                }
                
                for (int i = 0; i < turbineReader.RowCount; i++)
                {
                    Console.WriteLine("");
                    Console.Write("{0}              ", (i+1));
                    for (int j = 0; j < mastReader.RowCount; j++)
                    {
                        Console.Write("{0:0.00}             ", MastWeights[i, j]);
                    }
                }

                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("************** Optimised Uncertainty, Energy and Sigma*********************************************");
                Console.WriteLine("");
                Console.WriteLine("Energy: {0}", FinalEnergy);
                Console.WriteLine("Sigma: {0}", FinalSigma);
                Console.WriteLine("");
                Console.WriteLine("Uncertainty: {0}", FinalUncertainty);
                Console.WriteLine("");
                Console.WriteLine("****************************************************************************************************");

                int GetIndex(int[] ids, int id)
                {
                    for (int i = 0; i < ids.Length; i++)
                    {
                        if (ids[i] == id) return i;
                    }
                    return -1;

                }
            }, stackSize);
            th.Start();
            th.Join();
        }

    }
}






