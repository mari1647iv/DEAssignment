package sample;

import java.net.URL;
import java.util.ResourceBundle;
import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.chart.LineChart;
import javafx.scene.control.Label;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Separator;
import javafx.scene.layout.Pane;
import javafx.scene.chart.XYChart;
import javafx.scene.control.TextField;
import javafx.scene.paint.Color;

public class Controller implements Initializable {

    @FXML
    private LineChart<Number, Number> EMChart, IEMChart, RKMChart;
    @FXML
    private Label Numerator, Denominator, ErrorLabel;
    @FXML
    private Pane ErrorGraphData;
    @FXML
    private RadioButton ErrorButton, IVPButton;
    @FXML
    private Separator Sep1, Sep2;

    @FXML
    private TextField x0s, y0s, Xs, Ns, N2s;

    @Override
    public void initialize(URL url, ResourceBundle rb) {
        buildGraph();
    }

    public void buildGraph() {
        Func f = new Func();

        EMChart.getData().clear();
        IEMChart.getData().clear();
        RKMChart.getData().clear();
        double x0 = Double.parseDouble(x0s.getText()), y0 = Double.parseDouble(y0s.getText()), X = Double.parseDouble(Xs.getText());
        int N = Integer.parseInt(Ns.getText()), N2 = Integer.parseInt(N2s.getText());

        if (x0 > X) {
            ErrorLabel.setText("Wrong values of x0 and X! x0 must be less or equal X.");
            ErrorLabel.setTextFill(Color.RED);
            Numerator.setText("");
            Denominator.setText("");
        } else {
            if (x0 == 0) {
                ErrorLabel.setText("Wrong values of x0! This IVP cannot be solved.");
                ErrorLabel.setTextFill(Color.RED);
                Numerator.setText("");
                Denominator.setText("");
            } else {
                double c = (y0 * x0 + 2) / (y0 * Math.pow(x0, 5) - 2 * Math.pow(x0, 4));
                if ((x0 == 1) && (y0 == 0)) {
                    Numerator.setText("2*x^4 - 2");
                    Denominator.setText("x^5 + x");
                } else {
                    if (c < 0) {
                        if (c == -0.5) {
                            Numerator.setText("x^4 - 2");
                        } else {
                            Numerator.setText(-2 * c + "*x^4 - 2");
                        }
                        if (c == -1) {
                            Denominator.setText("x^5 + x");
                        } else {
                            Denominator.setText(-c + "*x^5 + x");
                        }
                    }
                    if (c == 0) {
                        Numerator.setText("-2");
                        Denominator.setText("x");
                    } else {
                        if (c == 0.5) {
                            Numerator.setText("x^4 + 2");
                        } else {
                            Numerator.setText(2 * c + "*x^4 + 2");
                        }
                        if (c == 1) {
                            Denominator.setText("x^5 - x");
                        } else {
                            Denominator.setText(c + "*x^5 - x");
                        }
                    }
                }
                if (N < 1) {
                        ErrorLabel.setText("Wrong value of N! Graph cannot be constructed.");
                        ErrorLabel.setTextFill(Color.RED);
                } else {
                    XYChart.Series<Number, Number> EM = new XYChart.Series<>();
                    XYChart.Series<Number, Number> IEM = new XYChart.Series<>();
                    XYChart.Series<Number, Number> RKM = new XYChart.Series<>();

                    if (IVPButton.isSelected()) {
                        XYChart.Series<Number, Number> EMTotalError = new XYChart.Series<>();
                        XYChart.Series<Number, Number> EMLocalError = new XYChart.Series<>();

                        XYChart.Series<Number, Number> IEMTotalError = new XYChart.Series<>();
                        XYChart.Series<Number, Number> IEMLocalError = new XYChart.Series<>();

                        XYChart.Series<Number, Number> RKMTotalError = new XYChart.Series<>();
                        XYChart.Series<Number, Number> RKMLocalError = new XYChart.Series<>();

                        XYChart.Series<Number, Number> ES = new XYChart.Series<>();
                        XYChart.Series<Number, Number> ESCopy1 = new XYChart.Series<>();
                        XYChart.Series<Number, Number> ESCopy2 = new XYChart.Series<>();

                        EM.setName("Euler Method Graph");
                        EMTotalError.setName("Euler Method Total Errors Graph");
                        EMLocalError.setName("Euler Method Local Errors Graph");

                        IEM.setName("Improved Euler Method Graph");
                        IEMTotalError.setName("Improved Euler Method Total Errors Graph");
                        IEMLocalError.setName("Improved Euler Method Local Errors Graph");

                        RKM.setName("Runge-Kutta Method Graph");
                        RKMTotalError.setName("Runge-Kutta Method Total Errors Graph");
                        RKMLocalError.setName("Runge-Kutta Method Local Errors Graph");

                        ES.setName("Exact Solution oF IVP");
                        ESCopy1.setName("Exact Solution oF IVP");
                        ESCopy2.setName("Exact Solution oF IVP");


                        double[] x = new double[N], y = new double[N], ete = new double[N], ele = new double[N], yi = new double[N], ite = new double[N], ile = new double[N], k = new double[4], yrk = new double[N], rkte = new double[N], rkle = new double[N], ye = new double[N];


                        x[0] = x0;
                        y[0] = y0;
                        yi[0] = y0;
                        yrk[0] = y0;
                        ye[0] = y0;

                        ete[0] = 0;
                        ite[0] = 0;
                        rkte[0] = 0;

                        ele[0] = 0;
                        ile[0] = 0;
                        rkle[0] = 0;

                        double h = (X - x0) / (N - 1);


                        int i0 = N;
                        for (int i = 1; i < N; i++) {
                            x[i] = x[i - 1] + h;
                            if ((x[i - 1] < 0) && (x[i] >= 0)) {
                                i0 = i;
                                break;
                            }
                            ye[i] = (2 * c * Math.pow(x[i], 4) + 2) / (c * Math.pow(x[i], 5) - x[i]);

                            y[i] = y[i - 1] + h * f.getResult(x[i - 1], y[i - 1]);
                            ete[i] = Math.abs(y[i] - ye[i]);
                            ele[i] = Math.abs(ete[i] - ete[i - 1]);

                            yi[i] = yi[i - 1] + (h / 2) * (f.getResult(x[i-1], yi[i-1]) + f.getResult(x[i], yi[i-1] + h * f.getResult(x[i-1], yi[i-1])));
                            ite[i] = Math.abs(yi[i] - ye[i]);
                            ile[i] = Math.abs(ite[i] - ite[i - 1]);

                            k[0] = f.getResult(x[i - 1], yrk[i - 1]);
                            k[1] = f.getResult(x[i - 1] + (h / 2), yrk[i - 1] + (k[0] * h / 2));
                            k[2] = f.getResult(x[i - 1] + (h / 2), yrk[i - 1] + (k[1] * h / 2));
                            k[3] = f.getResult(x[i - 1] + h, yrk[i - 1] + (h * k[2]));
                            yrk[i] = yrk[i - 1] + (h / 6) * (k[0] + 2 * k[1] + 2 * k[2] + k[3]);
                            rkte[i] = Math.abs(yrk[i] - ye[i]);
                            rkle[i] = Math.abs(rkte[i] - rkte[i - 1]);
                        }


                        for (int i = 0; i < i0; i++) {
                            EM.getData().add(new XYChart.Data<>(x[i], y[i]));
                            EMTotalError.getData().add(new XYChart.Data<>(x[i], ete[i]));
                            EMLocalError.getData().add(new XYChart.Data<>(x[i], ele[i]));

                            IEM.getData().add(new XYChart.Data<>(x[i], yi[i]));
                            IEMTotalError.getData().add(new XYChart.Data<>(x[i], ite[i]));
                            IEMLocalError.getData().add(new XYChart.Data<>(x[i], ile[i]));

                            RKM.getData().add(new XYChart.Data<>(x[i], yrk[i]));
                            RKMTotalError.getData().add(new XYChart.Data<>(x[i], rkte[i]));
                            RKMLocalError.getData().add(new XYChart.Data<>(x[i], rkle[i]));

                            ES.getData().add(new XYChart.Data<>(x[i], ye[i]));
                            ESCopy1.getData().add(new XYChart.Data<>(x[i], ye[i]));
                            ESCopy2.getData().add(new XYChart.Data<>(x[i], ye[i]));
                        }


                        EMChart.getData().setAll(EM);
                        EMChart.getData().addAll(ES, EMTotalError, EMLocalError);

                        IEMChart.getData().setAll(IEM);
                        IEMChart.getData().addAll(ESCopy2, IEMTotalError, IEMLocalError);

                        RKMChart.getData().setAll(RKM);
                        RKMChart.getData().addAll(ESCopy1, RKMTotalError, RKMLocalError);


                        if (i0 != N) {
                            ErrorLabel.setText("Graph was pertially constructed! y[0] does not exist. Graph cannot be fully constructed.");
                            ErrorLabel.setTextFill(Color.web("#FFDC00"));
                        } else {
                            ErrorLabel.setText("Graph was successfully constructed!");
                            ErrorLabel.setTextFill(Color.LIME);
                        }
                    }
                    if(ErrorButton.isSelected()) {

                        if (N < 1) {
                            ErrorLabel.setText("Wrong value of N1! Graph cannot be constructed.");
                            ErrorLabel.setTextFill(Color.RED);
                            EMChart.getData().clear();
                            IEMChart.getData().clear();
                            RKMChart.getData().clear();
                        } else {
                            if (N > N2) {
                                ErrorLabel.setText("Wrong values of N1 and N2! Graph cannot be constructed.");
                                ErrorLabel.setTextFill(Color.RED);
                                EMChart.getData().clear();
                                IEMChart.getData().clear();
                                RKMChart.getData().clear();
                            } else {
                                EM.setName("Total Error of Approximation in Euler Method");
                                IEM.setName("Total Error of Approximation in Improved Euler Method");
                                RKM.setName("Total Error of Approximation in Runge-Kutta Method");
                                double x, y, yi, yrk, yeX = (2 * c * Math.pow(X, 4) + 2) / (c * Math.pow(X, 5) - X);
                                double[] k = new double[4];

                                for (int i = 0; i < N2 - N + 1; i++) {
                                    double h = (X - x0) / (N + i - 1);

                                    x = x0;
                                    y = y0;
                                    yi = y0;
                                    yrk = y0;

                                    for (int j = 1; j < N + i; j++) {
                                        if ((x < 0) && (x + h >= 0)) {
                                            break;
                                        }

                                        k[0] = f.getResult(x, yrk);
                                        k[1] = f.getResult(x + (h / 2), yrk + (k[0] * h / 2));
                                        k[2] = f.getResult(x + (h / 2), yrk + (k[1] * h / 2));
                                        k[3] = f.getResult(x + h, yrk + (h * k[2]));

                                        y = y + h * f.getResult(x, y);
                                        x = x + h;
                                        yi = yi + (h / 2) * (f.getResult(x - h, yi) + f.getResult(x, yi + h * f.getResult(x - h, yi)));
                                        yrk = yrk + (h / 6) * (k[0] + 2 * k[1] + 2 * k[2] + k[3]);
                                    }

                                    EM.getData().add(new XYChart.Data<>(N + i, Math.abs(y - yeX)));
                                    IEM.getData().add(new XYChart.Data<>(N + i, Math.abs(yi - yeX)));
                                    RKM.getData().add(new XYChart.Data<>(N + i, Math.abs(yrk - yeX)));
                                }

                                IEMChart.getData().setAll(EM);
                                IEMChart.getData().addAll(IEM, RKM);

                                ErrorLabel.setText("Graph was successfully constructed!");
                                ErrorLabel.setTextFill(Color.LIME);
                            }
                        }
                    }
                }
            }
        }
    }

    public void toInitial() {
        x0s.setText("1");
        y0s.setText("0");
        Xs.setText("7");
        Ns.setText("25");
        N2s.setText("40");

        buildGraph();
    }

    public void ShowPane(){
        ErrorGraphData.setVisible(true);
        EMChart.setVisible(false);
        RKMChart.setVisible(false);
        Sep1.setVisible(false);
        Sep2.setVisible(false);
        IEMChart.setTitle("Errors of Approximation");

        toInitial();
    }
    public void HidePane(){
        ErrorGraphData.setVisible(false);
        EMChart.setVisible(true);
        RKMChart.setVisible(true);
        Sep1.setVisible(true);
        Sep2.setVisible(true);
        IEMChart.setTitle("Improved Euler Method");

        toInitial();
    }
}

class Func{

    public double getResult(double x, double y){
        return (4 / (x * x)) - (y / x) - (y * y);
    }
}

abstract class NumericalMethod{
    double x0, y0, X;
    int N;

    public XYChart.Series<Number, Number> calculateSeries;
}

class EulerMethod extends NumericalMethod{

    public EulerMethod(double X0, double Y0 , double x, int n) {
        x0 = X0;
        y0 = Y0;
        X = x;
        N = n;
    }

    //public XYChart.Series<Number, Number> calculateSeries(){

  //  }
}

class ImprovedEulerMethod extends NumericalMethod{

    public ImprovedEulerMethod(double X0, double Y0 , double x, int n) {
        x0 = X0;
        y0 = Y0;
        X = x;
        N = n;
    }

    //public XYChart.Series<Number, Number> calculateSeries(){

    //  }
}

class RungeKuttaMethod extends NumericalMethod{

    public RungeKuttaMethod(double X0, double Y0 , double x, int n) {
        x0 = X0;
        y0 = Y0;
        X = x;
        N = n;
    }

    //public XYChart.Series<Number, Number> calculateSeries(){

    //  }
}

class ExactSolution{
    double x0, y0, X, c;
    int N;

    public ExactSolution(double X0, double Y0 , double x, int n) {
        x0 = X0;
        y0 = Y0;
        X = x;
        N = n;
        c = (y0 * x0 + 2) / (y0 * Math.pow(x0, 5) - 2 * Math.pow(x0, 4));;
    }

    public double getY(double x){
        return (2 * c * Math.pow(x, 4) + 2) / (c * Math.pow(x, 5) - x);
    }

    //public XYChart.Series<Number, Number> calculateSeries(){

    //  }
}
