package sample;

import javafx.fxml.FXML;
import javafx.fxml.Initializable;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Label;
import javafx.scene.control.RadioButton;
import javafx.scene.control.Separator;
import javafx.scene.control.TextField;
import javafx.scene.layout.Pane;
import javafx.scene.paint.Color;

import java.net.URL;
import java.util.ResourceBundle;

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
                    EulerMethod EM;
                    ImprovedEulerMethod IEM;
                    RungeKuttaMethod RKM;
                    if (IVPButton.isSelected()) {
                        ExactSolution ES = new ExactSolution(x0, y0, X, N);

                        EM = new EulerMethod(x0, y0, X, N);
                        IEM = new ImprovedEulerMethod(x0, y0, X, N);
                        RKM = new RungeKuttaMethod(x0, y0, X, N);

                        EMChart.getData().setAll(EM.Y);
                        EMChart.getData().addAll(ES.calculateSeries(), EM.TotalError, EM.LocalError);

                        IEMChart.getData().setAll(IEM.Y);
                        IEMChart.getData().addAll(ES.calculateSeries(), IEM.TotalError, IEM.LocalError);

                        RKMChart.getData().setAll(RKM.Y);
                        RKMChart.getData().addAll(ES.calculateSeries(), RKM.TotalError, RKM.LocalError);


                        if ((x0 < 0) && (X >= 0)) {
                            ErrorLabel.setText("Graph was pertially constructed! y[0] does not exist. Graph cannot be fully constructed.");
                            ErrorLabel.setTextFill(Color.web("#FFDC00"));
                        } else {
                            ErrorLabel.setText("Graph was successfully constructed!");
                            ErrorLabel.setTextFill(Color.LIME);
                        }
                    }
                    if (ErrorButton.isSelected()) {

                        if (N < 1) {
                            ErrorLabel.setText("Wrong value of N1! Graph cannot be constructed.");
                            ErrorLabel.setTextFill(Color.RED);
                        } else {
                            if (N > N2) {
                                ErrorLabel.setText("Wrong values of N1 and N2! Graph cannot be constructed.");
                                ErrorLabel.setTextFill(Color.RED);
                            } else {
                                XYChart.Series<Number, Number> EMTotalError = new XYChart.Series<>();
                                XYChart.Series<Number, Number> IEMTotalError = new XYChart.Series<>();
                                XYChart.Series<Number, Number> RKMTotalError = new XYChart.Series<>();

                                EMTotalError.setName("Total Error of Approximation in Euler Method");
                                IEMTotalError.setName("Total Error of Approximation in Improved Euler Method");
                                RKMTotalError.setName("Total Error of Approximation in Runge-Kutta Method");

                                for (int i = 0; i < N2 - N + 1; i++) {
                                    EM = new EulerMethod(x0, y0, X, N + i);
                                    IEM = new ImprovedEulerMethod(x0, y0, X, N + i);
                                    RKM = new RungeKuttaMethod(x0, y0, X, N + i);

                                    EMTotalError.getData().add(new XYChart.Data<>(N + i, EM.TotalErrorInX));
                                    IEMTotalError.getData().add(new XYChart.Data<>(N + i, IEM.TotalErrorInX));
                                    RKMTotalError.getData().add(new XYChart.Data<>(N + i, RKM.TotalErrorInX));
                                }

                                IEMChart.getData().setAll(EMTotalError);
                                IEMChart.getData().addAll(IEMTotalError, RKMTotalError);

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

    public void ShowPane() {
        ErrorGraphData.setVisible(true);
        EMChart.setVisible(false);
        RKMChart.setVisible(false);
        Sep1.setVisible(false);
        Sep2.setVisible(false);
        IEMChart.setTitle("Errors of Approximation");

        toInitial();
    }

    public void HidePane() {
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
    double x0, y0, X, TotalErrorInX;
    int N;
    XYChart.Series<Number, Number> Y = new XYChart.Series<>();
    XYChart.Series<Number, Number> TotalError = new XYChart.Series<>();
    XYChart.Series<Number, Number> LocalError = new XYChart.Series<>();

    abstract public void calculateSeries();
}

class EulerMethod extends NumericalMethod{

    public EulerMethod(double X0, double Y0 , double x, int n) {
        x0 = X0;
        y0 = Y0;
        X = x;
        N = n;

        Y.setName("Euler Method Graph");
        TotalError.setName("Euler Method Total Errors Graph");
        LocalError.setName("Euler Method Local Errors Graph");

        calculateSeries();
    }

    public void calculateSeries(){
        Func f = new Func();
        ExactSolution es = new ExactSolution(x0, y0, X, N);
        double[] x = new double[N], y = new double[N], te = new double[N], le = new double[N];

        x[0] = x0;
        y[0] = y0;
        te[0] = 0;
        le[0] = 0;
        double h = (X - x0) / (N - 1);

        int i0 = N;
        for (int i = 1; i < N; i++) {
            x[i] = x[i - 1] + h;
            if ((x[i - 1] < 0) && (x[i] >= 0)) {
                i0 = i;
                break;
            }

            y[i] = y[i - 1] + h * f.getResult(x[i - 1], y[i - 1]);
            te[i] = Math.abs(y[i] - es.getY(x[i]));
            le[i] = Math.abs(te[i] - te[i - 1]);
        }
        TotalErrorInX = te[i0-1];

        for (int i = 0; i < i0; i++) {
            Y.getData().add(new XYChart.Data<>(x[i], y[i]));
            TotalError.getData().add(new XYChart.Data<>(x[i], te[i]));
            LocalError.getData().add(new XYChart.Data<>(x[i], le[i]));
        }
    }
}

class ImprovedEulerMethod extends NumericalMethod{

    public ImprovedEulerMethod(double X0, double Y0 , double x, int n) {
        x0 = X0;
        y0 = Y0;
        X = x;
        N = n;

        Y.setName("Improved Euler Method Graph");
        TotalError.setName("Improved Euler Method Total Errors Graph");
        LocalError.setName("Improved Euler Method Local Errors Graph");

        calculateSeries();
    }

    public void calculateSeries(){
        Func f = new Func();
        ExactSolution es = new ExactSolution(x0, y0, X, N);
        double[] x = new double[N], y = new double[N], te = new double[N], le = new double[N];

        x[0] = x0;
        y[0] = y0;
        te[0] = 0;
        le[0] = 0;
        double h = (X - x0) / (N - 1);

        int i0 = N;
        for (int i = 1; i < N; i++) {
            x[i] = x[i - 1] + h;
            if ((x[i - 1] < 0) && (x[i] >= 0)) {
                i0 = i;
                break;
            }

            y[i] = y[i - 1] + (h / 2) * (f.getResult(x[i-1], y[i-1]) + f.getResult(x[i], y[i-1] + h * f.getResult(x[i-1], y[i-1])));
            te[i] = Math.abs(y[i] - es.getY(x[i]));
            le[i] = Math.abs(te[i] - te[i - 1]);
        }

        TotalErrorInX = te[i0-1];

        for (int i = 0; i < i0; i++) {
            Y.getData().add(new XYChart.Data<>(x[i], y[i]));
            TotalError.getData().add(new XYChart.Data<>(x[i], te[i]));
            LocalError.getData().add(new XYChart.Data<>(x[i], le[i]));
        }
    }
}

class RungeKuttaMethod extends NumericalMethod{

    public RungeKuttaMethod(double X0, double Y0 , double x, int n) {
        x0 = X0;
        y0 = Y0;
        X = x;
        N = n;

        Y.setName("Runge-Kutta Method Graph");
        TotalError.setName("Runge-Kutta Method Total Errors Graph");
        LocalError.setName("Runge-Kutta Method Local Errors Graph");

        calculateSeries();
    }

    public void calculateSeries(){
        Func f = new Func();
        ExactSolution es = new ExactSolution(x0, y0, X, N);
        double[] x = new double[N], y = new double[N], te = new double[N], le = new double[N], k = new double[4];

        x[0] = x0;
        y[0] = y0;
        te[0] = 0;
        le[0] = 0;
        double h = (X - x0) / (N - 1);

        int i0 = N;
        for (int i = 1; i < N; i++) {
            x[i] = x[i - 1] + h;
            if ((x[i - 1] < 0) && (x[i] >= 0)) {
                i0 = i;
                break;
            }

            k[0] = f.getResult(x[i - 1], y[i - 1]);
            k[1] = f.getResult(x[i - 1] + (h / 2), y[i - 1] + (k[0] * h / 2));
            k[2] = f.getResult(x[i - 1] + (h / 2), y[i - 1] + (k[1] * h / 2));
            k[3] = f.getResult(x[i - 1] + h, y[i - 1] + (h * k[2]));
            y[i] = y[i - 1] + (h / 6) * (k[0] + 2 * k[1] + 2 * k[2] + k[3]);
            te[i] = Math.abs(y[i] - es.getY(x[i]));
            le[i] = Math.abs(te[i] - te[i - 1]);
        }

        TotalErrorInX = te[i0-1];

        for (int i = 0; i < i0; i++) {
            Y.getData().add(new XYChart.Data<>(x[i], y[i]));
            TotalError.getData().add(new XYChart.Data<>(x[i], te[i]));
            LocalError.getData().add(new XYChart.Data<>(x[i], le[i]));
        }
    }
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

    public XYChart.Series<Number, Number> calculateSeries(){
        XYChart.Series<Number, Number> Y = new XYChart.Series<>();
        Y.setName("Exact Solution of IVP");
        double x;

        double h = (X - x0) / (N - 1);
        x = x0;
        Y.getData().add(new XYChart.Data<>(x, y0));

        int i0 = N;
        for (int i = 1; i < N; i++) {
            if ((x < 0) && (x+h >= 0)) {
                i0 = i;
                break;
            }
            x = x + h;
            Y.getData().add(new XYChart.Data<>(x, getY(x)));
        }
        return Y;
    }
}
